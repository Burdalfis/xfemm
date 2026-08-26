#include "cudss_solver.h"

#include <cudss.h>
#include <cuda_runtime.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <queue>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace xfemm_benchmark {
namespace {

void cudaCheck(cudaError_t status, const char *expression)
{
    if (status != cudaSuccess)
        throw std::runtime_error(std::string(expression) + ": " +
                                 cudaGetErrorString(status));
}

void cudssCheck(cudssStatus_t status, const char *expression)
{
    if (status != CUDSS_STATUS_SUCCESS)
        throw std::runtime_error(std::string(expression) +
                                 ": cuDSS status " +
                                 std::to_string(static_cast<int>(status)));
}

#define CUDA_CHECK(expr) cudaCheck((expr), #expr)
#define CUDSS_CHECK(expr) cudssCheck((expr), #expr)

double elapsedMs(const std::chrono::steady_clock::time_point &start)
{
    return std::chrono::duration<double, std::milli>(
               std::chrono::steady_clock::now() - start).count();
}

template <typename Function>
double timeOnStream(cudaStream_t stream, Function &&function)
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

template <typename T>
class DeviceAllocation {
public:
    DeviceAllocation() = default;
    explicit DeviceAllocation(std::size_t count) : count_(count)
    {
        CUDA_CHECK(cudaMalloc(&pointer_, count_ * sizeof(T)));
    }
    DeviceAllocation(const DeviceAllocation &) = delete;
    DeviceAllocation &operator=(const DeviceAllocation &) = delete;
    ~DeviceAllocation() { if (pointer_) cudaFree(pointer_); }
    T *get() const { return pointer_; }
    std::size_t bytes() const { return count_ * sizeof(T); }
private:
    T *pointer_ = nullptr;
    std::size_t count_ = 0;
};

struct LowerCsr {
    std::vector<std::int32_t> rows;
    std::vector<std::int32_t> columns;
    std::vector<std::size_t> fullValueIndices;
};

LowerCsr makeLowerCsr(const LinearSystem &system)
{
    LowerCsr result;
    result.rows.resize(system.dimension() + 1);
    for (std::size_t row = 0; row < system.dimension(); ++row) {
        result.rows[row] = static_cast<std::int32_t>(result.columns.size());
        for (std::uint64_t j = system.rowOffsets[row];
             j < system.rowOffsets[row + 1]; ++j) {
            const std::int32_t column = system.columnIndices[j];
            if (column <= static_cast<std::int32_t>(row)) {
                result.columns.push_back(column);
                result.fullValueIndices.push_back(static_cast<std::size_t>(j));
            }
        }
    }
    result.rows.back() = static_cast<std::int32_t>(result.columns.size());
    return result;
}

std::vector<double> gatherLowerValues(const LinearSystem &system,
                                      const LowerCsr &lower)
{
    std::vector<double> result(lower.fullValueIndices.size());
    std::transform(lower.fullValueIndices.begin(), lower.fullValueIndices.end(),
                   result.begin(), [&](std::size_t index) {
                       return system.values[index];
                   });
    return result;
}

std::vector<std::uint64_t> undirectedOffDiagonalEdges(const LinearSystem &system)
{
    std::vector<std::uint64_t> edges;
    edges.reserve((system.nonzeros() - system.dimension()) / 2);
    for (std::size_t row = 0; row < system.dimension(); ++row) {
        for (std::uint64_t j = system.rowOffsets[row];
             j < system.rowOffsets[row + 1]; ++j) {
            const auto column = static_cast<std::size_t>(system.columnIndices[j]);
            if (row < column)
                edges.push_back((static_cast<std::uint64_t>(row) << 32) |
                                static_cast<std::uint64_t>(column));
        }
    }
    std::sort(edges.begin(), edges.end());
    return edges;
}

struct UniversalAgeTopology {
    std::vector<std::uint64_t> rows;
    std::vector<std::int32_t> columns;
    std::size_t changedEdges = 0;
    std::size_t components = 0;
    std::size_t firstRingNodes = 0;
    std::size_t secondRingNodes = 0;
};

UniversalAgeTopology buildUniversalAgeTopology(
    const LinearSystem &reference, const LinearSystem &shifted,
    const std::vector<LinearSystem> &systems, bool completeRevolutionSuperset)
{
    if (reference.dimension() != shifted.dimension())
        throw std::invalid_argument("AGE topology inputs have different dimensions");
    const auto firstEdges = undirectedOffDiagonalEdges(reference);
    const auto shiftedEdges = undirectedOffDiagonalEdges(shifted);
    std::vector<std::uint64_t> changed;
    std::set_symmetric_difference(firstEdges.begin(), firstEdges.end(),
                                  shiftedEdges.begin(), shiftedEdges.end(),
                                  std::back_inserter(changed));
    if (changed.empty())
        throw std::invalid_argument("AGE reference positions have identical topology");

    std::vector<std::vector<std::int32_t>> adjacency(reference.dimension());
    for (std::uint64_t edge : changed) {
        const auto a = static_cast<std::int32_t>(edge >> 32);
        const auto b = static_cast<std::int32_t>(edge & 0xffffffffu);
        adjacency[static_cast<std::size_t>(a)].push_back(b);
        adjacency[static_cast<std::size_t>(b)].push_back(a);
    }
    std::vector<signed char> color(reference.dimension(), -1);
    std::vector<std::int32_t> firstRing, secondRing;
    std::size_t components = 0;
    for (std::size_t start = 0; start < adjacency.size(); ++start) {
        if (adjacency[start].empty() || color[start] >= 0) continue;
        ++components;
        color[start] = 0;
        std::queue<std::int32_t> pending;
        pending.push(static_cast<std::int32_t>(start));
        while (!pending.empty()) {
            const auto node = pending.front();
            pending.pop();
            for (const auto neighbor : adjacency[static_cast<std::size_t>(node)]) {
                if (color[static_cast<std::size_t>(neighbor)] < 0) {
                    color[static_cast<std::size_t>(neighbor)] =
                        static_cast<signed char>(1 - color[static_cast<std::size_t>(node)]);
                    pending.push(neighbor);
                } else if (color[static_cast<std::size_t>(neighbor)] ==
                           color[static_cast<std::size_t>(node)]) {
                    throw std::runtime_error(
                        "changed AGE graph is not bipartite; cannot infer rings");
                }
            }
        }
    }
    for (std::size_t node = 0; node < color.size(); ++node) {
        if (color[node] == 0)
            firstRing.push_back(static_cast<std::int32_t>(node));
        else if (color[node] == 1)
            secondRing.push_back(static_cast<std::int32_t>(node));
    }
    if (firstRing.empty() || firstRing.size() != secondRing.size())
        throw std::runtime_error("could not infer two equal AGE node rings");
    if (completeRevolutionSuperset && components != 1)
        throw std::runtime_error(
            "AGE change graph is disconnected; choose reference positions "
            "farther apart to orient both complete-revolution rings");

    std::vector<std::vector<std::int32_t>> rowColumns(reference.dimension());
    for (const auto &system : systems) {
        if (system.dimension() != reference.dimension())
            throw std::invalid_argument("AGE systems have different dimensions");
        for (std::size_t row = 0; row < system.dimension(); ++row) {
            auto &target = rowColumns[row];
            for (std::uint64_t j = system.rowOffsets[row];
                 j < system.rowOffsets[row + 1]; ++j)
                target.push_back(system.columnIndices[j]);
        }
    }
    if (completeRevolutionSuperset) {
        // Across one complete relative revolution, every node on one uniformly
        // shifted ring can occupy every position relative to the other ring.
        // Retaining the complete cross product is therefore the exact structural
        // superset of all angle-specific cross-ring AGE couplings.
        for (const auto row : firstRing)
            rowColumns[static_cast<std::size_t>(row)].insert(
                rowColumns[static_cast<std::size_t>(row)].end(), secondRing.begin(),
                secondRing.end());
        for (const auto row : secondRing)
            rowColumns[static_cast<std::size_t>(row)].insert(
                rowColumns[static_cast<std::size_t>(row)].end(), firstRing.begin(),
                firstRing.end());
    }

    UniversalAgeTopology result;
    result.rows.resize(reference.dimension() + 1);
    for (std::size_t row = 0; row < rowColumns.size(); ++row) {
        auto &columns = rowColumns[row];
        std::sort(columns.begin(), columns.end());
        columns.erase(std::unique(columns.begin(), columns.end()), columns.end());
        result.rows[row] = result.columns.size();
        result.columns.insert(result.columns.end(), columns.begin(), columns.end());
    }
    result.rows.back() = result.columns.size();
    result.changedEdges = changed.size();
    result.components = components;
    result.firstRingNodes = firstRing.size();
    result.secondRingNodes = secondRing.size();
    return result;
}

LinearSystem padToTopology(const LinearSystem &source,
                           const UniversalAgeTopology &topology)
{
    LinearSystem result = source;
    result.sourcePath += "#universal-age";
    result.rowOffsets = topology.rows;
    result.columnIndices = topology.columns;
    result.values.assign(topology.columns.size(), 0.);
    for (std::size_t row = 0; row < source.dimension(); ++row) {
        std::size_t target = static_cast<std::size_t>(topology.rows[row]);
        const std::size_t targetEnd = static_cast<std::size_t>(topology.rows[row + 1]);
        for (std::uint64_t j = source.rowOffsets[row];
             j < source.rowOffsets[row + 1]; ++j) {
            const auto column = source.columnIndices[j];
            while (target < targetEnd && topology.columns[target] < column)
                ++target;
            if (target == targetEnd || topology.columns[target] != column)
                throw std::runtime_error("universal AGE topology omitted an input entry");
            result.values[target] = source.values[j];
        }
    }
    return result;
}

class CudssContext {
public:
    explicit CudssContext(const LinearSystem &system)
        : topology_(topologyHash(system)), lower_(makeLowerCsr(system)),
          dRows_(lower_.rows.size()), dColumns_(lower_.columns.size()),
          dValues_(lower_.columns.size()), dRhs_(system.dimension()),
          dSolution_(system.dimension())
    {
        const auto setupStarted = std::chrono::steady_clock::now();
        CUDSS_CHECK(cudssCreate(&handle_));
        CUDA_CHECK(cudaStreamCreateWithFlags(&stream_, cudaStreamNonBlocking));
        CUDSS_CHECK(cudssSetStream(handle_, stream_));
        CUDSS_CHECK(cudssConfigCreate(&config_));
        const char *deterministic = std::getenv("XFEMM_CUDSS_DETERMINISTIC");
        if (deterministic && deterministic[0] != '\0' && deterministic[0] != '0') {
            const int enabled = 1;
            CUDSS_CHECK(cudssConfigSet(config_, CUDSS_CONFIG_DETERMINISTIC_MODE,
                                       &enabled, sizeof(enabled)));
        }
        CUDSS_CHECK(cudssDataCreate(handle_, &data_));
        const std::int64_t n = static_cast<std::int64_t>(system.dimension());
        const std::int64_t nnz = static_cast<std::int64_t>(lower_.columns.size());
        CUDSS_CHECK(cudssMatrixCreateCsr(
            &matrix_, n, n, nnz, dRows_.get(), nullptr,
            dColumns_.get(), dValues_.get(), CUDSS_R_32I, CUDSS_R_32I,
            CUDSS_R_64F, CUDSS_MTYPE_SPD, CUDSS_MVIEW_LOWER,
            CUDSS_BASE_ZERO));
        CUDSS_CHECK(cudssMatrixCreateDn(&rhs_, n, 1, n, dRhs_.get(),
                                        CUDSS_R_64F, CUDSS_LAYOUT_COL_MAJOR));
        CUDSS_CHECK(cudssMatrixCreateDn(&solution_, n, 1, n,
                                        dSolution_.get(), CUDSS_R_64F,
                                        CUDSS_LAYOUT_COL_MAJOR));
        setupMs_ = elapsedMs(setupStarted);

        const auto values = gatherLowerValues(system, lower_);
        uploadMs_ = timeOnStream(stream_, [&] {
            CUDA_CHECK(cudaMemcpyAsync(dRows_.get(), lower_.rows.data(), dRows_.bytes(),
                                       cudaMemcpyHostToDevice, stream_));
            CUDA_CHECK(cudaMemcpyAsync(dColumns_.get(), lower_.columns.data(),
                                       dColumns_.bytes(), cudaMemcpyHostToDevice,
                                       stream_));
            CUDA_CHECK(cudaMemcpyAsync(dValues_.get(), values.data(), dValues_.bytes(),
                                       cudaMemcpyHostToDevice, stream_));
            CUDA_CHECK(cudaMemcpyAsync(dRhs_.get(), system.rhs.data(), dRhs_.bytes(),
                                       cudaMemcpyHostToDevice, stream_));
            CUDA_CHECK(cudaMemsetAsync(dSolution_.get(), 0, dSolution_.bytes(),
                                       stream_));
        });
        analysisMs_ = timeOnStream(stream_, [&] {
            CUDSS_CHECK(cudssExecute(handle_, CUDSS_PHASE_ANALYSIS, config_, data_,
                                     matrix_, solution_, rhs_));
        });
        std::size_t written = 0;
        const cudssStatus_t memoryStatus = cudssDataGet(
            handle_, data_, CUDSS_DATA_MEMORY_ESTIMATES,
            memoryEstimates_.data(), memoryEstimates_.size() * sizeof(std::int64_t),
            &written);
        if (memoryStatus != CUDSS_STATUS_SUCCESS)
            memoryEstimates_.fill(0);
    }

    CudssContext(const CudssContext &) = delete;
    ~CudssContext()
    {
        if (solution_) cudssMatrixDestroy(solution_);
        if (rhs_) cudssMatrixDestroy(rhs_);
        if (matrix_) cudssMatrixDestroy(matrix_);
        if (data_) cudssDataDestroy(handle_, data_);
        if (config_) cudssConfigDestroy(config_);
        if (handle_) cudssDestroy(handle_);
        if (stream_) cudaStreamDestroy(stream_);
    }

    BenchmarkResult solve(const LinearSystem &system, bool first,
                          bool validate = true,
                          std::vector<double> *solutionOut = nullptr)
    {
        if (topologyHash(system) != topology_)
            throw std::invalid_argument("cuDSS value update changed CSR topology");
        if (!first) {
            const auto values = gatherLowerValues(system, lower_);
            uploadMs_ = timeOnStream(stream_, [&] {
                CUDA_CHECK(cudaMemcpyAsync(dValues_.get(), values.data(),
                                           dValues_.bytes(), cudaMemcpyHostToDevice,
                                           stream_));
                CUDA_CHECK(cudaMemcpyAsync(dRhs_.get(), system.rhs.data(),
                                           dRhs_.bytes(), cudaMemcpyHostToDevice,
                                           stream_));
            });
        }

        const int phase = first ? CUDSS_PHASE_FACTORIZATION
                                : CUDSS_PHASE_REFACTORIZATION;
        const double factorMs = timeOnStream(stream_, [&] {
            CUDSS_CHECK(cudssExecute(handle_, phase, config_, data_, matrix_,
                                     solution_, rhs_));
        });
        const double solveMs = timeOnStream(stream_, [&] {
            CUDSS_CHECK(cudssExecute(handle_, CUDSS_PHASE_SOLVE, config_, data_,
                                     matrix_, solution_, rhs_));
        });

        std::vector<double> solution;
        if (validate || solutionOut) {
            solution.resize(system.dimension());
            CUDA_CHECK(cudaMemcpyAsync(solution.data(), dSolution_.get(),
                                       dSolution_.bytes(), cudaMemcpyDeviceToHost,
                                       stream_));
            CUDA_CHECK(cudaStreamSynchronize(stream_));
            if (solutionOut) *solutionOut = solution;
        }
        std::int64_t factorNnz = 0;
        std::size_t written = 0;
        if (cudssDataGet(handle_, data_, CUDSS_DATA_LU_NNZ, &factorNnz,
                         sizeof(factorNnz), &written) != CUDSS_STATUS_SUCCESS)
            factorNnz = 0;

        BenchmarkResult result;
        result.system = system.sourcePath;
        result.method = first ? "cuDSS" : "cuDSS-refactor";
        result.preconditioner = "direct-SPD";
        result.dimension = system.dimension();
        result.nonzeros = system.nonzeros();
        result.setupMs = first ? setupMs_ : 0.;
        result.uploadMs = uploadMs_;
        result.analysisMs = first ? analysisMs_ : 0.;
        result.factorizationMs = factorMs;
        result.solveMs = solveMs;
        result.factorNonzeros = static_cast<std::uint64_t>(std::max<std::int64_t>(0, factorNnz));
        result.permanentMemoryBytes = static_cast<std::uint64_t>(
            std::max<std::int64_t>(0, memoryEstimates_[0]));
        result.memoryBytes = static_cast<std::uint64_t>(
            std::max<std::int64_t>(0, memoryEstimates_[1]));
        if (validate) {
            result.relativeResidual = relativeResidual(system, solution);
            result.solutionError = relativeSolutionError(solution,
                                                         system.cpuSolution);
            result.solutionHash = solutionHash(solution);
            if (referenceSolution_.empty() ||
                referenceSolveIndex_ != system.solveIndex) {
                referenceSolution_ = solution;
                referenceSolveIndex_ = system.solveIndex;
            } else
                result.repeatabilityError =
                    relativeSolutionError(solution, referenceSolution_);
        }
        return result;
    }

private:
    std::uint64_t topology_ = 0;
    LowerCsr lower_;
    DeviceAllocation<std::int32_t> dRows_, dColumns_;
    DeviceAllocation<double> dValues_, dRhs_, dSolution_;
    cudssHandle_t handle_ = nullptr;
    cudaStream_t stream_ = nullptr;
    cudssConfig_t config_ = nullptr;
    cudssData_t data_ = nullptr;
    cudssMatrix_t matrix_ = nullptr;
    cudssMatrix_t rhs_ = nullptr;
    cudssMatrix_t solution_ = nullptr;
    std::array<std::int64_t, 16> memoryEstimates_{};
    double setupMs_ = 0.;
    double uploadMs_ = 0.;
    double analysisMs_ = 0.;
    std::vector<double> referenceSolution_;
    std::uint64_t referenceSolveIndex_ = ~std::uint64_t{0};
};

} // namespace

std::vector<BenchmarkResult>
runCudssSequence(const std::vector<LinearSystem> &systems, int repetitions)
{
    std::vector<BenchmarkResult> results;
    if (systems.empty()) return results;
    CudssContext context(systems.front());
    results.push_back(context.solve(systems.front(), true));
    for (int repeat = 1; repeat < repetitions; ++repeat) {
        auto result = context.solve(systems.front(), true);
        result.setupMs = 0.;
        result.uploadMs = 0.;
        result.analysisMs = 0.;
        results.push_back(result);
    }
    const std::uint64_t topology = topologyHash(systems.front());
    for (std::size_t i = 1; i < systems.size(); ++i) {
        if (topologyHash(systems[i]) != topology) break;
        results.push_back(context.solve(systems[i], false));
    }
    return results;
}

std::vector<CudssBatchResult>
runCudssConcurrency(const LinearSystem &system,
                    const std::vector<int> &batchSizes)
{
    std::vector<CudssBatchResult> output;
    if (batchSizes.empty()) return output;
    const int largest = *std::max_element(batchSizes.begin(), batchSizes.end());
    std::vector<std::unique_ptr<CudssContext>> contexts;
    contexts.reserve(static_cast<std::size_t>(largest));
    for (int i = 0; i < largest; ++i)
        contexts.emplace_back(new CudssContext(system));

    for (int batch : batchSizes) {
        std::atomic<bool> start{false};
        std::vector<BenchmarkResult> results(static_cast<std::size_t>(batch));
        std::vector<std::thread> workers;
        workers.reserve(static_cast<std::size_t>(batch));
        for (int i = 0; i < batch; ++i) {
            workers.emplace_back([&, i] {
                while (!start.load(std::memory_order_acquire))
                    std::this_thread::yield();
                results[static_cast<std::size_t>(i)] =
                    contexts[static_cast<std::size_t>(i)]->solve(system, true,
                                                                 false);
            });
        }
        const auto wallStarted = std::chrono::steady_clock::now();
        start.store(true, std::memory_order_release);
        for (auto &worker : workers) worker.join();

        CudssBatchResult result;
        result.batch = batch;
        result.wallMs = elapsedMs(wallStarted);
        for (const auto &item : results) {
            result.averageFactorizationMs += item.factorizationMs;
            result.averageSolveMs += item.solveMs;
            result.memoryBytes += item.memoryBytes;
        }
        result.averageFactorizationMs /= batch;
        result.averageSolveMs /= batch;
        output.push_back(result);
    }
    return output;
}

CudssAgeTopologyResult
runCudssAgeTopologyExperiment(const LinearSystem &reference,
                              const LinearSystem &shifted,
                              const std::vector<LinearSystem> &systems,
                              bool completeRevolutionSuperset)
{
    if (systems.empty())
        throw std::invalid_argument("AGE topology experiment needs systems");
    const auto buildStarted = std::chrono::steady_clock::now();
    const auto topology = buildUniversalAgeTopology(reference, shifted, systems,
                                                     completeRevolutionSuperset);
    std::vector<LinearSystem> universalSystems;
    universalSystems.reserve(systems.size());
    for (const auto &system : systems)
        universalSystems.push_back(padToTopology(system, topology));

    CudssAgeTopologyResult result;
    result.universalBuildMs = elapsedMs(buildStarted);
    result.changedUndirectedEdges = topology.changedEdges;
    result.changedGraphComponents = topology.components;
    result.firstRingNodes = topology.firstRingNodes;
    result.secondRingNodes = topology.secondRingNodes;
    result.angleSpecificNonzeros = systems.front().nonzeros();
    result.universalNonzeros = topology.columns.size();
    const auto universalLower = makeLowerCsr(universalSystems.front());
    result.universalLowerNonzeros = universalLower.columns.size();
    result.universalUpdateBytes =
        universalLower.columns.size() * sizeof(double) +
        systems.front().dimension() * sizeof(double);

    const auto specificStarted = std::chrono::steady_clock::now();
    std::unordered_map<std::uint64_t, std::unique_ptr<CudssContext>>
        specificContexts;
    std::uint64_t specificMaximumTemporaryBytes = 0;
    for (const auto &system : systems) {
        const auto hash = topologyHash(system);
        auto context = specificContexts.find(hash);
        const bool first = context == specificContexts.end();
        if (first) {
            context = specificContexts.emplace(
                hash, std::unique_ptr<CudssContext>(new CudssContext(system))).first;
        }
        result.angleSpecific.push_back(context->second->solve(system, first));
        if (first) {
            const auto &measurement = result.angleSpecific.back();
            result.angleSpecificCacheResidentBytes +=
                measurement.permanentMemoryBytes;
            const auto temporaryBytes = measurement.memoryBytes >
                                        measurement.permanentMemoryBytes
                                      ? measurement.memoryBytes -
                                        measurement.permanentMemoryBytes
                                      : 0;
            specificMaximumTemporaryBytes = std::max(
                specificMaximumTemporaryBytes, temporaryBytes);
        }
    }
    result.angleSpecificCachePeakBytes =
        result.angleSpecificCacheResidentBytes + specificMaximumTemporaryBytes;
    result.angleSpecificWallMs = elapsedMs(specificStarted);

    const auto universalStarted = std::chrono::steady_clock::now();
    CudssContext context(universalSystems.front());
    result.universalSolutions.resize(universalSystems.size());
    result.universal.push_back(context.solve(universalSystems.front(), true, true,
                                              &result.universalSolutions.front()));
    result.universalResidentBytes = result.universal.front().permanentMemoryBytes;
    result.universalPeakBytes = result.universal.front().memoryBytes;
    for (std::size_t i = 1; i < universalSystems.size(); ++i)
        result.universal.push_back(context.solve(universalSystems[i], false, true,
                                                  &result.universalSolutions[i]));
    result.universalWallMs = elapsedMs(universalStarted);
    return result;
}

} // namespace xfemm_benchmark
