#include "xfemm_system.h"

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusparse.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace xb = xfemm_benchmark;

namespace {

void check(cudaError_t status, const char *expression)
{
    if (status != cudaSuccess)
        throw std::runtime_error(std::string(expression) + ": " +
                                 cudaGetErrorString(status));
}
void check(cublasStatus_t status, const char *expression)
{
    if (status != CUBLAS_STATUS_SUCCESS)
        throw std::runtime_error(std::string(expression) + ": cuBLAS error");
}
void check(cusparseStatus_t status, const char *expression)
{
    if (status != CUSPARSE_STATUS_SUCCESS)
        throw std::runtime_error(std::string(expression) + ": cuSPARSE error " +
                                 std::to_string(static_cast<int>(status)));
}
#define CUDA_CHECK(expr) check((expr), #expr)
#define CUBLAS_CHECK(expr) check((expr), #expr)
#define CUSPARSE_CHECK(expr) check((expr), #expr)

template <typename T>
class DeviceBuffer {
public:
    DeviceBuffer() = default;
    explicit DeviceBuffer(std::size_t count) : count_(count)
    {
        CUDA_CHECK(cudaMalloc(&data_, count_ * sizeof(T)));
    }
    DeviceBuffer(const DeviceBuffer &) = delete;
    DeviceBuffer &operator=(const DeviceBuffer &) = delete;
    DeviceBuffer(DeviceBuffer &&other) noexcept
        : data_(other.data_), count_(other.count_)
    {
        other.data_ = nullptr;
        other.count_ = 0;
    }
    DeviceBuffer &operator=(DeviceBuffer &&other) noexcept
    {
        if (this != &other) {
            if (data_) cudaFree(data_);
            data_ = other.data_;
            count_ = other.count_;
            other.data_ = nullptr;
            other.count_ = 0;
        }
        return *this;
    }
    ~DeviceBuffer() { if (data_) cudaFree(data_); }
    T *get() const { return data_; }
    std::size_t bytes() const { return count_ * sizeof(T); }
private:
    T *data_ = nullptr;
    std::size_t count_ = 0;
};

template <typename Function>
double timeGpu(cudaStream_t stream, Function &&function)
{
    cudaEvent_t start = nullptr, stop = nullptr;
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

__global__ void prepareSsorRows(int n, const std::int32_t *rows,
                                const std::int32_t *columns,
                                const double *matrixValues, double lambda,
                                double *lowerValues, double *diagonal)
{
    const int row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row >= n) return;
    double d = 0.;
    for (int j = rows[row]; j < rows[row + 1]; ++j) {
        const int column = columns[j];
        if (column < row)
            lowerValues[j] = lambda * matrixValues[j];
        else if (column == row) {
            d = matrixValues[j];
            lowerValues[j] = d;
        } else {
            lowerValues[j] = 0.;
        }
    }
    diagonal[row] = d;
}

__global__ void multiplyDiagonal(int n, const double *diagonal, double *values)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) values[i] *= diagonal[i];
}

class SsorPcg {
public:
    explicit SsorPcg(const xb::LinearSystem &system)
        : topology_(xb::topologyHash(system)), n_(static_cast<int>(system.dimension())),
          nnz_(static_cast<int>(system.nonzeros())), rows_(system.dimension() + 1),
          columns_(system.nonzeros()), matrixValues_(system.nonzeros()),
          lowerValues_(system.nonzeros()), diagonal_(system.dimension()),
          b_(system.dimension()), x_(system.dimension()), r_(system.dimension()),
          p_(system.dimension()), z_(system.dimension()), ap_(system.dimension()),
          work_(system.dimension())
    {
        std::vector<std::int32_t> compactRows(system.rowOffsets.size());
        std::transform(system.rowOffsets.begin(), system.rowOffsets.end(),
                       compactRows.begin(), [](std::uint64_t value) {
                           return static_cast<std::int32_t>(value);
                       });
        CUDA_CHECK(cudaStreamCreateWithFlags(&stream_, cudaStreamNonBlocking));
        CUBLAS_CHECK(cublasCreate(&blas_));
        CUSPARSE_CHECK(cusparseCreate(&sparse_));
        CUBLAS_CHECK(cublasSetStream(blas_, stream_));
        CUSPARSE_CHECK(cusparseSetStream(sparse_, stream_));
        structureUploadMs_ = timeGpu(stream_, [&] {
            CUDA_CHECK(cudaMemcpyAsync(rows_.get(), compactRows.data(), rows_.bytes(),
                                       cudaMemcpyHostToDevice, stream_));
            CUDA_CHECK(cudaMemcpyAsync(columns_.get(), system.columnIndices.data(),
                                       columns_.bytes(), cudaMemcpyHostToDevice,
                                       stream_));
        });
        CUSPARSE_CHECK(cusparseCreateCsr(
            &matrix_, n_, n_, nnz_, rows_.get(), columns_.get(),
            matrixValues_.get(), CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
            CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F));
        CUSPARSE_CHECK(cusparseCreateCsr(
            &lower_, n_, n_, nnz_, rows_.get(), columns_.get(), lowerValues_.get(),
            CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO,
            CUDA_R_64F));
        cusparseFillMode_t fillMode = CUSPARSE_FILL_MODE_LOWER;
        cusparseDiagType_t diagonalType = CUSPARSE_DIAG_TYPE_NON_UNIT;
        CUSPARSE_CHECK(cusparseSpMatSetAttribute(
            lower_, CUSPARSE_SPMAT_FILL_MODE, &fillMode, sizeof(fillMode)));
        CUSPARSE_CHECK(cusparseSpMatSetAttribute(
            lower_, CUSPARSE_SPMAT_DIAG_TYPE, &diagonalType,
            sizeof(diagonalType)));
        CUSPARSE_CHECK(cusparseCreateDnVec(&vecP_, n_, p_.get(), CUDA_R_64F));
        CUSPARSE_CHECK(cusparseCreateDnVec(&vecAp_, n_, ap_.get(), CUDA_R_64F));
        CUSPARSE_CHECK(cusparseCreateDnVec(&forwardInput_, n_, r_.get(), CUDA_R_64F));
        CUSPARSE_CHECK(cusparseCreateDnVec(&forwardOutput_, n_, work_.get(), CUDA_R_64F));
        CUSPARSE_CHECK(cusparseCreateDnVec(&backwardInput_, n_, work_.get(), CUDA_R_64F));
        CUSPARSE_CHECK(cusparseCreateDnVec(&backwardOutput_, n_, z_.get(), CUDA_R_64F));
        CUSPARSE_CHECK(cusparseSpSV_createDescr(&forwardDescription_));
        CUSPARSE_CHECK(cusparseSpSV_createDescr(&backwardDescription_));

        std::size_t spmvBytes = 0, forwardBytes = 0, backwardBytes = 0;
        CUSPARSE_CHECK(cusparseSpMV_bufferSize(
            sparse_, CUSPARSE_OPERATION_NON_TRANSPOSE, &one_, matrix_, vecP_,
            &zero_, vecAp_, CUDA_R_64F, CUSPARSE_SPMV_CSR_ALG1, &spmvBytes));
        CUSPARSE_CHECK(cusparseSpSV_bufferSize(
            sparse_, CUSPARSE_OPERATION_NON_TRANSPOSE, &ssorScale_, lower_,
            forwardInput_, forwardOutput_, CUDA_R_64F, CUSPARSE_SPSV_ALG_DEFAULT,
            forwardDescription_, &forwardBytes));
        CUSPARSE_CHECK(cusparseSpSV_bufferSize(
            sparse_, CUSPARSE_OPERATION_TRANSPOSE, &one_, lower_, backwardInput_,
            backwardOutput_, CUDA_R_64F, CUSPARSE_SPSV_ALG_DEFAULT,
            backwardDescription_, &backwardBytes));
        spmvBuffer_ = DeviceBuffer<unsigned char>(spmvBytes);
        // SpSV analysis stores persistent information in its external buffer.
        // The forward and transpose descriptors therefore need independent
        // buffers; sharing one lets the second analysis overwrite the first.
        forwardBuffer_ = DeviceBuffer<unsigned char>(forwardBytes);
        backwardBuffer_ = DeviceBuffer<unsigned char>(backwardBytes);
        upload(system, true);
    }

    ~SsorPcg()
    {
        if (backwardDescription_) cusparseSpSV_destroyDescr(backwardDescription_);
        if (forwardDescription_) cusparseSpSV_destroyDescr(forwardDescription_);
        if (backwardOutput_) cusparseDestroyDnVec(backwardOutput_);
        if (backwardInput_) cusparseDestroyDnVec(backwardInput_);
        if (forwardOutput_) cusparseDestroyDnVec(forwardOutput_);
        if (forwardInput_) cusparseDestroyDnVec(forwardInput_);
        if (vecAp_) cusparseDestroyDnVec(vecAp_);
        if (vecP_) cusparseDestroyDnVec(vecP_);
        if (lower_) cusparseDestroySpMat(lower_);
        if (matrix_) cusparseDestroySpMat(matrix_);
        if (sparse_) cusparseDestroy(sparse_);
        if (blas_) cublasDestroy(blas_);
        if (stream_) cudaStreamDestroy(stream_);
    }

    void upload(const xb::LinearSystem &system, bool uploadInitial)
    {
        if (xb::topologyHash(system) != topology_)
            throw std::invalid_argument("SSOR update changed topology");
        const double valueUploadMs = timeGpu(stream_, [&] {
            CUDA_CHECK(cudaMemcpyAsync(matrixValues_.get(), system.values.data(),
                                       matrixValues_.bytes(), cudaMemcpyHostToDevice,
                                       stream_));
            CUDA_CHECK(cudaMemcpyAsync(b_.get(), system.rhs.data(), b_.bytes(),
                                       cudaMemcpyHostToDevice, stream_));
            if (uploadInitial)
                CUDA_CHECK(cudaMemcpyAsync(x_.get(), system.initialSolution.data(),
                                           x_.bytes(), cudaMemcpyHostToDevice,
                                           stream_));
        });
        uploadMs_ = valueUploadMs + (structureUploadPending_ ? structureUploadMs_
                                                              : 0.);
        preparationMs_ = timeGpu(stream_, [&] {
            const int blocks = (n_ + 255) / 256;
            prepareSsorRows<<<blocks, 256, 0, stream_>>>(
                n_, rows_.get(), columns_.get(), matrixValues_.get(), lambda_,
                lowerValues_.get(), diagonal_.get());
            CUDA_CHECK(cudaGetLastError());
            if (analyzed_) {
                CUSPARSE_CHECK(cusparseSpSV_updateMatrix(
                    sparse_, forwardDescription_, lowerValues_.get(),
                    CUSPARSE_SPSV_UPDATE_GENERAL));
                CUSPARSE_CHECK(cusparseSpSV_updateMatrix(
                    sparse_, backwardDescription_, lowerValues_.get(),
                    CUSPARSE_SPSV_UPDATE_GENERAL));
            }
        });
        if (!analyzed_) {
            analysisMs_ = timeGpu(stream_, [&] {
                CUSPARSE_CHECK(cusparseSpMV_preprocess(
                    sparse_, CUSPARSE_OPERATION_NON_TRANSPOSE, &one_, matrix_,
                    vecP_, &zero_, vecAp_, CUDA_R_64F, CUSPARSE_SPMV_CSR_ALG1,
                    spmvBuffer_.get()));
                CUSPARSE_CHECK(cusparseSpSV_analysis(
                    sparse_, CUSPARSE_OPERATION_NON_TRANSPOSE, &ssorScale_, lower_,
                    forwardInput_, forwardOutput_, CUDA_R_64F,
                    CUSPARSE_SPSV_ALG_DEFAULT, forwardDescription_,
                    forwardBuffer_.get()));
                CUSPARSE_CHECK(cusparseSpSV_analysis(
                    sparse_, CUSPARSE_OPERATION_TRANSPOSE, &one_, lower_,
                    backwardInput_, backwardOutput_, CUDA_R_64F,
                    CUSPARSE_SPSV_ALG_DEFAULT, backwardDescription_,
                    backwardBuffer_.get()));
            });
            analyzed_ = true;
        } else {
            analysisMs_ = 0.;
        }
    }

    void solve(const xb::LinearSystem &system, bool resetInitial,
               const char *record)
    {
        if (resetInitial)
            CUDA_CHECK(cudaMemcpyAsync(x_.get(), system.initialSolution.data(),
                                       x_.bytes(), cudaMemcpyHostToDevice, stream_));
        iterations_ = 0;
        stoppingMetric_ = 0.;
        const double solveMs = timeGpu(stream_, [&] { pcg(system.tolerance); });
        std::vector<double> solution(system.dimension());
        CUDA_CHECK(cudaMemcpyAsync(solution.data(), x_.get(), x_.bytes(),
                                   cudaMemcpyDeviceToHost, stream_));
        CUDA_CHECK(cudaStreamSynchronize(stream_));
        std::cout << record << ',' << system.sourcePath
                  << ",cuSPARSE-PCG,SpSV-SSOR," << system.dimension() << ','
                  << system.nonzeros() << ',' << std::fixed << std::setprecision(6)
                  << uploadMs_ << ',' << analysisMs_ << ',' << preparationMs_ << ','
                  << solveMs << ',' << iterations_ << ',' << memoryBytes() << ','
                  << std::scientific << std::setprecision(9) << stoppingMetric_ << ','
                  << xb::relativeResidual(system, solution) << ','
                  << xb::relativeSolutionError(solution, system.cpuSolution)
                  << ",0," 
                  << std::hex << xb::solutionHash(solution) << std::dec << '\n';
        structureUploadPending_ = false;
    }

private:
    void applyPreconditioner(double *input, double *output)
    {
        CUSPARSE_CHECK(cusparseDnVecSetValues(forwardInput_, input));
        CUSPARSE_CHECK(cusparseDnVecSetValues(backwardOutput_, output));
        CUSPARSE_CHECK(cusparseSpSV_solve(
            sparse_, CUSPARSE_OPERATION_NON_TRANSPOSE, &ssorScale_, lower_,
            forwardInput_, forwardOutput_, CUDA_R_64F,
            CUSPARSE_SPSV_ALG_DEFAULT, forwardDescription_));
        const int blocks = (n_ + 255) / 256;
        multiplyDiagonal<<<blocks, 256, 0, stream_>>>(n_, diagonal_.get(),
                                                      work_.get());
        CUDA_CHECK(cudaGetLastError());
        CUSPARSE_CHECK(cusparseSpSV_solve(
            sparse_, CUSPARSE_OPERATION_TRANSPOSE, &one_, lower_, backwardInput_,
            backwardOutput_, CUDA_R_64F, CUSPARSE_SPSV_ALG_DEFAULT,
            backwardDescription_));
    }

    void spmv()
    {
        CUSPARSE_CHECK(cusparseSpMV(
            sparse_, CUSPARSE_OPERATION_NON_TRANSPOSE, &one_, matrix_, vecP_,
            &zero_, vecAp_, CUDA_R_64F, CUSPARSE_SPMV_CSR_ALG1,
            spmvBuffer_.get()));
    }

    void pcg(double tolerance)
    {
        applyPreconditioner(b_.get(), z_.get());
        double initialNorm = 0.;
        CUBLAS_CHECK(cublasDdot(blas_, n_, z_.get(), 1, b_.get(), 1,
                                &initialNorm));
        if (initialNorm == 0.) return;
        CUBLAS_CHECK(cublasDcopy(blas_, n_, b_.get(), 1, r_.get(), 1));
        CUSPARSE_CHECK(cusparseDnVecSetValues(vecP_, x_.get()));
        spmv();
        const double negativeOne = -1.;
        CUBLAS_CHECK(cublasDaxpy(blas_, n_, &negativeOne, ap_.get(), 1,
                                r_.get(), 1));
        applyPreconditioner(r_.get(), z_.get());
        CUBLAS_CHECK(cublasDcopy(blas_, n_, z_.get(), 1, p_.get(), 1));
        CUSPARSE_CHECK(cusparseDnVecSetValues(vecP_, p_.get()));
        double residual = 0.;
        CUBLAS_CHECK(cublasDdot(blas_, n_, z_.get(), 1, r_.get(), 1, &residual));
        while (iterations_ < 10000) {
            spmv();
            double pAp = 0.;
            CUBLAS_CHECK(cublasDdot(blas_, n_, p_.get(), 1, ap_.get(), 1, &pAp));
            const double alpha = residual / pAp;
            CUBLAS_CHECK(cublasDaxpy(blas_, n_, &alpha, p_.get(), 1, x_.get(), 1));
            const double negativeAlpha = -alpha;
            CUBLAS_CHECK(cublasDaxpy(blas_, n_, &negativeAlpha, ap_.get(), 1,
                                    r_.get(), 1));
            applyPreconditioner(r_.get(), z_.get());
            double nextResidual = 0.;
            CUBLAS_CHECK(cublasDdot(blas_, n_, z_.get(), 1, r_.get(), 1,
                                    &nextResidual));
            stoppingMetric_ = std::sqrt(std::abs(nextResidual / initialNorm));
            ++iterations_;
            if (iterations_ % 1000 == 0)
                std::cerr << "SpSV-SSOR iterations=" << iterations_
                          << " stopping_metric=" << stoppingMetric_ << '\n';
            if (!std::isfinite(stoppingMetric_) || stoppingMetric_ <= tolerance)
                break;
            const double beta = nextResidual / residual;
            residual = nextResidual;
            CUBLAS_CHECK(cublasDscal(blas_, n_, &beta, p_.get(), 1));
            CUBLAS_CHECK(cublasDaxpy(blas_, n_, &one_, z_.get(), 1, p_.get(), 1));
        }
    }

    std::uint64_t memoryBytes() const
    {
        return rows_.bytes() + columns_.bytes() + matrixValues_.bytes() +
               lowerValues_.bytes() + diagonal_.bytes() + b_.bytes() +
               x_.bytes() + r_.bytes() + p_.bytes() + z_.bytes() + ap_.bytes() +
               work_.bytes() + spmvBuffer_.bytes() + forwardBuffer_.bytes() +
               backwardBuffer_.bytes();
    }

    std::uint64_t topology_ = 0;
    int n_ = 0, nnz_ = 0;
    DeviceBuffer<std::int32_t> rows_, columns_;
    DeviceBuffer<double> matrixValues_, lowerValues_, diagonal_;
    DeviceBuffer<double> b_, x_, r_, p_, z_, ap_, work_;
    DeviceBuffer<unsigned char> spmvBuffer_, forwardBuffer_, backwardBuffer_;
    cudaStream_t stream_ = nullptr;
    cublasHandle_t blas_ = nullptr;
    cusparseHandle_t sparse_ = nullptr;
    cusparseSpMatDescr_t matrix_ = nullptr, lower_ = nullptr;
    cusparseDnVecDescr_t vecP_ = nullptr, vecAp_ = nullptr;
    cusparseDnVecDescr_t forwardInput_ = nullptr, forwardOutput_ = nullptr;
    cusparseDnVecDescr_t backwardInput_ = nullptr, backwardOutput_ = nullptr;
    cusparseSpSVDescr_t forwardDescription_ = nullptr;
    cusparseSpSVDescr_t backwardDescription_ = nullptr;
    bool analyzed_ = false;
    double uploadMs_ = 0., structureUploadMs_ = 0., analysisMs_ = 0.;
    double preparationMs_ = 0.;
    bool structureUploadPending_ = true;
    std::int64_t iterations_ = 0;
    double stoppingMetric_ = 0.;
    const double lambda_ = 1.5;
    const double ssorScale_ = lambda_ * (2. - lambda_);
    const double one_ = 1., zero_ = 0.;
};

} // namespace

int main(int argc, char **argv)
{
    if (argc < 2) {
        std::cerr << "usage: xfemm-cuda-ssor-benchmark SYSTEM [SYSTEM ...]\n";
        return 2;
    }
    try {
        std::vector<xb::LinearSystem> systems;
        for (int i = 1; i < argc; ++i) systems.push_back(xb::readLinearSystem(argv[i]));
        std::cout << "record,system,method,preconditioner,n,nnz,upload_ms,"
                     "structural_analysis_ms,value_preparation_ms,solve_ms,"
                     "iterations,memory_bytes,stopping_metric,relative_residual,"
                     "solution_error,repeatability_error,solution_hash\n";
        SsorPcg solver(systems.front());
        solver.solve(systems.front(), true, "single");
        const auto topology = xb::topologyHash(systems.front());
        for (std::size_t i = 1; i < systems.size(); ++i) {
            if (xb::topologyHash(systems[i]) != topology) break;
            solver.upload(systems[i], false);
            solver.solve(systems[i], false, "value-update");
        }
    } catch (const std::exception &error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
