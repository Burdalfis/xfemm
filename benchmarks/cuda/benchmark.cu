#include "benchmark_result.h"
#include "xfemm_system.h"

#ifdef XFEMM_HAVE_CUDSS
#include "cudss_solver.h"
#endif

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusparse.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace xb = xfemm_benchmark;

namespace {

const char *cudaErrorName(cudaError_t status) { return cudaGetErrorString(status); }
const char *cudaErrorName(cublasStatus_t) { return "cuBLAS error"; }
const char *cudaErrorName(cusparseStatus_t) { return "cuSPARSE error"; }

template <typename Status>
void check(Status status, Status success, const char *expression)
{
    if (status != success)
        throw std::runtime_error(std::string(expression) + ": " +
                                 cudaErrorName(status));
}

#define CUDA_CHECK(expr) check((expr), cudaSuccess, #expr)
#define CUBLAS_CHECK(expr) check((expr), CUBLAS_STATUS_SUCCESS, #expr)
#define CUSPARSE_CHECK(expr) check((expr), CUSPARSE_STATUS_SUCCESS, #expr)

template <typename T>
struct DeviceBuffer {
    T *pointer = nullptr;
    std::size_t count = 0;
    DeviceBuffer() = default;
    explicit DeviceBuffer(std::size_t size) : count(size)
    {
        CUDA_CHECK(cudaMalloc(&pointer, count * sizeof(T)));
    }
    DeviceBuffer(const DeviceBuffer &) = delete;
    DeviceBuffer &operator=(const DeviceBuffer &) = delete;
    DeviceBuffer(DeviceBuffer &&other) noexcept
        : pointer(other.pointer), count(other.count)
    {
        other.pointer = nullptr;
        other.count = 0;
    }
    DeviceBuffer &operator=(DeviceBuffer &&other) noexcept
    {
        if (this != &other) {
            if (pointer) cudaFree(pointer);
            pointer = other.pointer;
            count = other.count;
            other.pointer = nullptr;
            other.count = 0;
        }
        return *this;
    }
    ~DeviceBuffer() { if (pointer) cudaFree(pointer); }
    std::size_t bytes() const { return count * sizeof(T); }
};

double elapsedMilliseconds(const std::chrono::steady_clock::time_point &start)
{
    return std::chrono::duration<double, std::milli>(
               std::chrono::steady_clock::now() - start).count();
}

template <typename Function>
double timeGpu(cudaStream_t stream, Function &&function)
{
    cudaEvent_t start = nullptr;
    cudaEvent_t stop = nullptr;
    CUDA_CHECK(cudaEventCreate(&start));
    CUDA_CHECK(cudaEventCreate(&stop));
    CUDA_CHECK(cudaEventRecord(start, stream));
    function();
    CUDA_CHECK(cudaEventRecord(stop, stream));
    CUDA_CHECK(cudaEventSynchronize(stop));
    float milliseconds = 0.;
    CUDA_CHECK(cudaEventElapsedTime(&milliseconds, start, stop));
    cudaEventDestroy(stop);
    cudaEventDestroy(start);
    return milliseconds;
}

__global__ void extractDiagonal(int n, const std::int32_t *rowOffsets,
                                const std::int32_t *columns,
                                const double *values, double *diagonal)
{
    const int row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row >= n) return;
    double value = 0.;
    for (int j = rowOffsets[row]; j < rowOffsets[row + 1]; ++j)
        if (columns[j] == row) {
            value = values[j];
            break;
        }
    diagonal[row] = value;
}

__global__ void applyJacobi(int n, const double *input,
                            const double *diagonal, double *output)
{
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < n) output[index] = input[index] / diagonal[index];
}

struct DeviceStructure {
    DeviceBuffer<std::int32_t> rows;
    DeviceBuffer<std::int32_t> columns;
    double uploadMs = 0.;
    std::uint64_t hash = 0;

    explicit DeviceStructure(const xb::LinearSystem &system)
        : rows(system.dimension() + 1), columns(system.nonzeros()),
          hash(xb::topologyHash(system))
    {
        std::vector<std::int32_t> compactRows(system.rowOffsets.size());
        std::transform(system.rowOffsets.begin(), system.rowOffsets.end(),
                       compactRows.begin(), [](std::uint64_t value) {
                           return static_cast<std::int32_t>(value);
                       });
        uploadMs = timeGpu(nullptr, [&] {
            CUDA_CHECK(cudaMemcpyAsync(rows.pointer, compactRows.data(), rows.bytes(),
                                       cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpyAsync(columns.pointer,
                                       system.columnIndices.data(), columns.bytes(),
                                       cudaMemcpyHostToDevice));
        });
    }
};

class JacobiPcg {
public:
    JacobiPcg(std::shared_ptr<DeviceStructure> structure,
              const xb::LinearSystem &system)
        : structure_(std::move(structure)), n_(static_cast<int>(system.dimension())),
          nnz_(static_cast<int>(system.nonzeros())), values_(system.nonzeros()),
          b_(system.dimension()), x_(system.dimension()), r_(system.dimension()),
          p_(system.dimension()), z_(system.dimension()), ap_(system.dimension()),
          diagonal_(system.dimension())
    {
        const auto started = std::chrono::steady_clock::now();
        CUDA_CHECK(cudaStreamCreateWithFlags(&stream_, cudaStreamNonBlocking));
        CUBLAS_CHECK(cublasCreate(&blas_));
        CUSPARSE_CHECK(cusparseCreate(&sparse_));
        CUBLAS_CHECK(cublasSetStream(blas_, stream_));
        CUSPARSE_CHECK(cusparseSetStream(sparse_, stream_));
        CUSPARSE_CHECK(cusparseCreateCsr(
            &matrix_, n_, n_, nnz_, structure_->rows.pointer,
            structure_->columns.pointer, values_.pointer,
            CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO,
            CUDA_R_64F));
        CUSPARSE_CHECK(cusparseCreateDnVec(&vectorP_, n_, p_.pointer, CUDA_R_64F));
        CUSPARSE_CHECK(cusparseCreateDnVec(&vectorAp_, n_, ap_.pointer, CUDA_R_64F));
        CUSPARSE_CHECK(cusparseSpMV_bufferSize(
            sparse_, CUSPARSE_OPERATION_NON_TRANSPOSE, &one_, matrix_, vectorP_,
            &zero_, vectorAp_, CUDA_R_64F, CUSPARSE_SPMV_CSR_ALG1,
            &spmvBufferBytes_));
        spmvBuffer_ = DeviceBuffer<unsigned char>(spmvBufferBytes_);
        setupMs_ = elapsedMilliseconds(started);
        upload(system, true);
    }

    JacobiPcg(const JacobiPcg &) = delete;
    ~JacobiPcg()
    {
        if (vectorAp_) cusparseDestroyDnVec(vectorAp_);
        if (vectorP_) cusparseDestroyDnVec(vectorP_);
        if (matrix_) cusparseDestroySpMat(matrix_);
        if (sparse_) cusparseDestroy(sparse_);
        if (blas_) cublasDestroy(blas_);
        if (stream_) cudaStreamDestroy(stream_);
    }

    void upload(const xb::LinearSystem &system, bool useExportedInitial)
    {
        if (xb::topologyHash(system) != structure_->hash)
            throw std::invalid_argument("PCG value update changed CSR topology");
        uploadMs_ = timeGpu(stream_, [&] {
            CUDA_CHECK(cudaMemcpyAsync(values_.pointer, system.values.data(),
                                       values_.bytes(), cudaMemcpyHostToDevice,
                                       stream_));
            CUDA_CHECK(cudaMemcpyAsync(b_.pointer, system.rhs.data(), b_.bytes(),
                                       cudaMemcpyHostToDevice, stream_));
            if (useExportedInitial)
                CUDA_CHECK(cudaMemcpyAsync(x_.pointer,
                                           system.initialSolution.data(), x_.bytes(),
                                           cudaMemcpyHostToDevice, stream_));
        });
        preparationMs_ = timeGpu(stream_, [&] {
            const int blocks = (n_ + 255) / 256;
            extractDiagonal<<<blocks, 256, 0, stream_>>>(
                n_, structure_->rows.pointer, structure_->columns.pointer,
                values_.pointer, diagonal_.pointer);
            CUDA_CHECK(cudaGetLastError());
            CUSPARSE_CHECK(cusparseSpMV_preprocess(
                sparse_, CUSPARSE_OPERATION_NON_TRANSPOSE, &one_, matrix_, vectorP_,
                &zero_, vectorAp_, CUDA_R_64F, CUSPARSE_SPMV_CSR_ALG1,
                spmvBuffer_.pointer));
        });
    }

    xb::BenchmarkResult solve(const xb::LinearSystem &system,
                              bool resetToExportedInitial)
    {
        if (resetToExportedInitial)
            CUDA_CHECK(cudaMemcpyAsync(x_.pointer, system.initialSolution.data(),
                                       x_.bytes(), cudaMemcpyHostToDevice, stream_));
        iterations_ = 0;
        stoppingMetric_ = 0.;
        const double solveMs = timeGpu(stream_, [&] { solveInternal(system.tolerance); });

        std::vector<double> solution(system.dimension());
        CUDA_CHECK(cudaMemcpyAsync(solution.data(), x_.pointer, x_.bytes(),
                                   cudaMemcpyDeviceToHost, stream_));
        CUDA_CHECK(cudaStreamSynchronize(stream_));

        xb::BenchmarkResult result;
        result.system = system.sourcePath;
        result.method = "cuSPARSE-PCG";
        result.preconditioner = "Jacobi";
        result.dimension = system.dimension();
        result.nonzeros = system.nonzeros();
        result.setupMs = setupPending_ ? setupMs_ : 0.;
        result.uploadMs = uploadMs_ +
                          (structureUploadPending_ ? structure_->uploadMs : 0.);
        result.preparationMs = preparationMs_;
        result.solveMs = solveMs;
        result.iterations = iterations_;
        result.memoryBytes = memoryBytes();
        result.stoppingMetric = stoppingMetric_;
        result.relativeResidual = xb::relativeResidual(system, solution);
        result.solutionError =
            xb::relativeSolutionError(solution, system.cpuSolution);
        result.solutionHash = xb::solutionHash(solution);
        if (referenceSolution_.empty() || referenceSolveIndex_ != system.solveIndex) {
            referenceSolution_ = solution;
            referenceSolveIndex_ = system.solveIndex;
        } else
            result.repeatabilityError =
                xb::relativeSolutionError(solution, referenceSolution_);
        setupPending_ = false;
        structureUploadPending_ = false;
        uploadMs_ = 0.;
        preparationMs_ = 0.;
        return result;
    }

private:
    void spmv()
    {
        CUSPARSE_CHECK(cusparseSpMV(
            sparse_, CUSPARSE_OPERATION_NON_TRANSPOSE, &one_, matrix_, vectorP_,
            &zero_, vectorAp_, CUDA_R_64F, CUSPARSE_SPMV_CSR_ALG1,
            spmvBuffer_.pointer));
    }

    void solveInternal(double tolerance)
    {
        const int blocks = (n_ + 255) / 256;
        applyJacobi<<<blocks, 256, 0, stream_>>>(n_, b_.pointer,
                                                diagonal_.pointer, z_.pointer);
        CUDA_CHECK(cudaGetLastError());
        double initialNorm = 0.;
        CUBLAS_CHECK(cublasDdot(blas_, n_, z_.pointer, 1, b_.pointer, 1,
                                &initialNorm));
        if (initialNorm == 0.) return;

        CUBLAS_CHECK(cublasDcopy(blas_, n_, b_.pointer, 1, r_.pointer, 1));
        // Repoint the reusable SpMV descriptors for the initial residual.
        CUSPARSE_CHECK(cusparseDnVecSetValues(vectorP_, x_.pointer));
        CUSPARSE_CHECK(cusparseDnVecSetValues(vectorAp_, ap_.pointer));
        spmv();
        const double minusOne = -1.;
        CUBLAS_CHECK(cublasDaxpy(blas_, n_, &minusOne, ap_.pointer, 1,
                                r_.pointer, 1));
        applyJacobi<<<blocks, 256, 0, stream_>>>(n_, r_.pointer,
                                                diagonal_.pointer, z_.pointer);
        CUDA_CHECK(cudaGetLastError());
        CUBLAS_CHECK(cublasDcopy(blas_, n_, z_.pointer, 1, p_.pointer, 1));
        double residual = 0.;
        CUBLAS_CHECK(cublasDdot(blas_, n_, z_.pointer, 1, r_.pointer, 1,
                                &residual));

        CUSPARSE_CHECK(cusparseDnVecSetValues(vectorP_, p_.pointer));
        while (iterations_ < 100000) {
            spmv();
            double pAp = 0.;
            CUBLAS_CHECK(cublasDdot(blas_, n_, p_.pointer, 1, ap_.pointer, 1,
                                    &pAp));
            const double alpha = residual / pAp;
            CUBLAS_CHECK(cublasDaxpy(blas_, n_, &alpha, p_.pointer, 1,
                                    x_.pointer, 1));
            const double negativeAlpha = -alpha;
            CUBLAS_CHECK(cublasDaxpy(blas_, n_, &negativeAlpha, ap_.pointer, 1,
                                    r_.pointer, 1));
            applyJacobi<<<blocks, 256, 0, stream_>>>(
                n_, r_.pointer, diagonal_.pointer, z_.pointer);
            CUDA_CHECK(cudaGetLastError());
            double nextResidual = 0.;
            CUBLAS_CHECK(cublasDdot(blas_, n_, z_.pointer, 1, r_.pointer, 1,
                                    &nextResidual));
            stoppingMetric_ = std::sqrt(std::abs(nextResidual / initialNorm));
            ++iterations_;
            if (!std::isfinite(stoppingMetric_) || stoppingMetric_ <= tolerance)
                break;
            const double beta = nextResidual / residual;
            residual = nextResidual;
            CUBLAS_CHECK(cublasDscal(blas_, n_, &beta, p_.pointer, 1));
            CUBLAS_CHECK(cublasDaxpy(blas_, n_, &one_, z_.pointer, 1,
                                    p_.pointer, 1));
        }
    }

    std::uint64_t memoryBytes() const
    {
        return structure_->rows.bytes() + structure_->columns.bytes() +
               values_.bytes() + b_.bytes() + x_.bytes() + r_.bytes() +
               p_.bytes() + z_.bytes() + ap_.bytes() + diagonal_.bytes() +
               spmvBuffer_.bytes();
    }

    std::shared_ptr<DeviceStructure> structure_;
    int n_ = 0;
    int nnz_ = 0;
    cudaStream_t stream_ = nullptr;
    cublasHandle_t blas_ = nullptr;
    cusparseHandle_t sparse_ = nullptr;
    cusparseSpMatDescr_t matrix_ = nullptr;
    cusparseDnVecDescr_t vectorP_ = nullptr;
    cusparseDnVecDescr_t vectorAp_ = nullptr;
    DeviceBuffer<double> values_, b_, x_, r_, p_, z_, ap_, diagonal_;
    DeviceBuffer<unsigned char> spmvBuffer_;
    std::size_t spmvBufferBytes_ = 0;
    double setupMs_ = 0.;
    double uploadMs_ = 0.;
    double preparationMs_ = 0.;
    bool setupPending_ = true;
    bool structureUploadPending_ = true;
    std::int64_t iterations_ = 0;
    double stoppingMetric_ = 0.;
    std::vector<double> referenceSolution_;
    std::uint64_t referenceSolveIndex_ = ~std::uint64_t{0};
    const double one_ = 1.;
    const double zero_ = 0.;
};

void printHeader()
{
    std::cout << "record,system,method,preconditioner,n,nnz,setup_ms,upload_ms,"
                 "analysis_ms,preparation_ms,factorization_ms,solve_ms,iterations,"
                 "memory_bytes,factor_nnz,stopping_metric,relative_residual,"
                 "solution_error,repeatability_error,solution_hash\n";
}

void printResult(const char *record, const xb::BenchmarkResult &r)
{
    std::cout << record << ',' << r.system << ',' << r.method << ','
              << r.preconditioner << ',' << r.dimension << ',' << r.nonzeros
              << ',' << std::fixed << std::setprecision(6) << r.setupMs << ','
              << r.uploadMs << ',' << r.analysisMs << ',' << r.preparationMs
              << ',' << r.factorizationMs << ',' << r.solveMs << ','
              << r.iterations << ',' << r.memoryBytes << ',' << r.factorNonzeros
              << ',' << std::scientific << std::setprecision(9)
              << r.stoppingMetric << ',' << r.relativeResidual << ','
              << r.solutionError << ',' << r.repeatabilityError << ','
              << std::hex << r.solutionHash << std::dec << '\n';
}

std::vector<int> parseBatchSizes(const std::string &text)
{
    std::vector<int> result;
    std::stringstream stream(text);
    std::string item;
    while (std::getline(stream, item, ',')) {
        const int value = std::stoi(item);
        if (value <= 0 || value > 256) throw std::invalid_argument("invalid batch size");
        result.push_back(value);
    }
    return result;
}

void runConcurrency(const xb::LinearSystem &system,
                    const std::vector<int> &batchSizes)
{
    for (int batch : batchSizes) {
        auto structure = std::make_shared<DeviceStructure>(system);
        std::vector<std::unique_ptr<JacobiPcg>> solvers;
        for (int i = 0; i < batch; ++i)
            solvers.emplace_back(new JacobiPcg(structure, system));
        // Warm JIT/library caches once before the measured batch.
        solvers.front()->solve(system, true);

        std::atomic<bool> start{false};
        std::vector<xb::BenchmarkResult> results(static_cast<std::size_t>(batch));
        std::vector<std::thread> workers;
        for (int i = 0; i < batch; ++i)
            workers.emplace_back([&, i] {
                while (!start.load(std::memory_order_acquire))
                    std::this_thread::yield();
                results[static_cast<std::size_t>(i)] =
                    solvers[static_cast<std::size_t>(i)]->solve(system, true);
            });
        const auto wallStarted = std::chrono::steady_clock::now();
        start.store(true, std::memory_order_release);
        for (auto &worker : workers) worker.join();
        const double wallMs = elapsedMilliseconds(wallStarted);
        double sumLatency = 0.;
        std::uint64_t memory = structure->rows.bytes() + structure->columns.bytes();
        for (const auto &result : results) {
            sumLatency += result.solveMs;
            memory += result.memoryBytes - structure->rows.bytes() -
                      structure->columns.bytes();
        }
        std::cout << "batch," << system.sourcePath << ",cuSPARSE-PCG,Jacobi,"
                  << system.dimension() << ',' << system.nonzeros() << ','
                  << batch << ',' << std::fixed << std::setprecision(6)
                  << wallMs << ',' << (1000. * batch / wallMs) << ','
                  << 0. << ',' << (sumLatency / batch) << ',' << memory << '\n';
    }
}

} // namespace

int main(int argc, char **argv)
{
    try {
        std::vector<std::string> paths;
        std::vector<int> batchSizes;
        int repetitions = 3;
        bool runPcg = true;
        bool runCudss = true;
        for (int i = 1; i < argc; ++i) {
            const std::string argument = argv[i];
            if (argument.rfind("--repetitions=", 0) == 0)
                repetitions = std::stoi(argument.substr(14));
            else if (argument.rfind("--batch=", 0) == 0)
                batchSizes = parseBatchSizes(argument.substr(8));
            else if (argument == "--method=pcg") runCudss = false;
            else if (argument == "--method=cudss") runPcg = false;
            else if (!argument.empty() && argument[0] == '-')
                throw std::invalid_argument("unknown option: " + argument);
            else paths.push_back(argument);
        }
        if (paths.empty()) {
            std::cerr << "usage: xfemm-cuda-benchmark [--method=pcg|cudss] "
                         "[--repetitions=N] [--batch=1,2,...] SYSTEM...\n";
            return 2;
        }

        cudaDeviceProp properties{};
        CUDA_CHECK(cudaGetDeviceProperties(&properties, 0));
        std::cerr << "GPU name=" << properties.name << " compute="
                  << properties.major << '.' << properties.minor
                  << " memory_bytes=" << properties.totalGlobalMem << '\n';

        std::vector<xb::LinearSystem> systems;
        for (const auto &path : paths) systems.push_back(xb::readLinearSystem(path));
        printHeader();
        if (runPcg) {
            auto structure = std::make_shared<DeviceStructure>(systems.front());
            JacobiPcg solver(structure, systems.front());
            for (int repeat = 0; repeat < repetitions; ++repeat)
                printResult("single", solver.solve(systems.front(), true));

            // Same-topology values/RHS updates, retaining the preceding GPU x.
            for (std::size_t i = 1; i < systems.size(); ++i) {
                if (xb::topologyHash(systems[i]) != structure->hash) break;
                solver.upload(systems[i], false);
                printResult("value-update", solver.solve(systems[i], false));
            }
        }
#ifdef XFEMM_HAVE_CUDSS
        if (runCudss)
            for (const auto &result : xb::runCudssSequence(systems, repetitions))
                printResult("direct", result);
#else
        if (runCudss)
            std::cerr << "cuDSS was not available at build time\n";
#endif
        if (!batchSizes.empty()) {
            std::cout << "record,system,method,preconditioner,n,nnz,batch,wall_ms,"
                         "systems_per_second,average_factorization_ms,"
                         "average_solve_ms,memory_bytes\n";
            if (runPcg) runConcurrency(systems.front(), batchSizes);
#ifdef XFEMM_HAVE_CUDSS
            if (runCudss) {
                for (const auto &result :
                     xb::runCudssConcurrency(systems.front(), batchSizes)) {
                    std::cout << "batch," << systems.front().sourcePath
                              << ",cuDSS,direct-SPD,"
                              << systems.front().dimension() << ','
                              << systems.front().nonzeros() << ',' << result.batch
                              << ',' << std::fixed << std::setprecision(6)
                              << result.wallMs << ','
                              << (1000. * result.batch / result.wallMs) << ','
                              << result.averageFactorizationMs << ','
                              << result.averageSolveMs << ','
                              << result.memoryBytes << '\n';
                }
            }
#endif
        }
    } catch (const std::exception &error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
