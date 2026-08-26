#ifndef XFEMM_CUDSS_SOLVER_H
#define XFEMM_CUDSS_SOLVER_H

#include "benchmark_result.h"
#include "xfemm_system.h"

#include <vector>

namespace xfemm_benchmark {

std::vector<BenchmarkResult>
runCudssSequence(const std::vector<LinearSystem> &systems, int repetitions);

struct CudssBatchResult {
    int batch = 0;
    double wallMs = 0.;
    double averageFactorizationMs = 0.;
    double averageSolveMs = 0.;
    std::uint64_t memoryBytes = 0;
};

struct CudssAgeTopologyResult {
    std::size_t changedUndirectedEdges = 0;
    std::size_t changedGraphComponents = 0;
    std::size_t firstRingNodes = 0;
    std::size_t secondRingNodes = 0;
    std::size_t angleSpecificNonzeros = 0;
    std::size_t universalNonzeros = 0;
    std::size_t universalLowerNonzeros = 0;
    std::uint64_t universalUpdateBytes = 0;
    std::uint64_t angleSpecificCacheResidentBytes = 0;
    std::uint64_t angleSpecificCachePeakBytes = 0;
    std::uint64_t universalResidentBytes = 0;
    std::uint64_t universalPeakBytes = 0;
    double universalBuildMs = 0.;
    double angleSpecificWallMs = 0.;
    double universalWallMs = 0.;
    /** One persistent context is cached per distinct topology. */
    std::vector<BenchmarkResult> angleSpecific;
    std::vector<BenchmarkResult> universal;
    /** Captured to permit end-to-end validation through XFEMM's existing
     * postprocessor.  The benchmark executable only writes these on request. */
    std::vector<std::vector<double>> universalSolutions;
};

std::vector<CudssBatchResult>
runCudssConcurrency(const LinearSystem &system,
                    const std::vector<int> &batchSizes);

/** Build the complete cross-ring AGE superset inferred from two positions,
 * pad every supplied numeric system with explicit zeros, and compare a
 * persistent per-topology cuDSS cache against one retained union analysis.
 * Set completeRevolutionSuperset=false to union only supplied patterns. */
CudssAgeTopologyResult
runCudssAgeTopologyExperiment(const LinearSystem &reference,
                              const LinearSystem &shifted,
                              const std::vector<LinearSystem> &systems,
                              bool completeRevolutionSuperset = true);

} // namespace xfemm_benchmark

#endif
