#ifndef XFEMM_CUDA_BENCHMARK_RESULT_H
#define XFEMM_CUDA_BENCHMARK_RESULT_H

#include <cstddef>
#include <cstdint>
#include <string>

namespace xfemm_benchmark {

struct BenchmarkResult {
    std::string system;
    std::string method;
    std::string preconditioner;
    std::size_t dimension = 0;
    std::size_t nonzeros = 0;
    double setupMs = 0.;
    double uploadMs = 0.;
    double analysisMs = 0.;
    double preparationMs = 0.;
    double factorizationMs = 0.;
    double solveMs = 0.;
    std::int64_t iterations = 0;
    /** cuDSS permanent device allocation estimate after analysis. */
    std::uint64_t permanentMemoryBytes = 0;
    /** Peak device allocation estimate (or total bytes for iterative methods). */
    std::uint64_t memoryBytes = 0;
    std::uint64_t factorNonzeros = 0;
    std::uint64_t solutionHash = 0;
    double stoppingMetric = 0.;
    double relativeResidual = 0.;
    double solutionError = 0.;
    double repeatabilityError = 0.;
};

} // namespace xfemm_benchmark

#endif
