#ifndef XFEMM_BENCHMARK_SYSTEM_H
#define XFEMM_BENCHMARK_SYSTEM_H

#include <cstdint>
#include <string>
#include <vector>

namespace xfemm_benchmark {

struct LinearSystem {
    std::string sourcePath;
    std::uint64_t solveIndex = 0;
    std::uint32_t flags = 0;
    std::uint32_t cpuThreads = 0;
    double tolerance = 0.;
    double ssorRelaxation = 0.;
    std::int64_t cpuIterations = 0;
    double cpuPcgResidual = 0.;
    std::vector<std::uint64_t> rowOffsets;
    std::vector<std::int32_t> columnIndices;
    std::vector<double> values;
    std::vector<double> rhs;
    std::vector<double> initialSolution;
    std::vector<double> cpuSolution;

    std::size_t dimension() const { return rhs.size(); }
    std::size_t nonzeros() const { return values.size(); }
    bool warmStart() const { return (flags & 1u) != 0; }
};

LinearSystem readLinearSystem(const std::string &path);
std::uint64_t topologyHash(const LinearSystem &system);
std::uint64_t solutionHash(const std::vector<double> &solution);
double relativeResidual(const LinearSystem &system,
                        const std::vector<double> &solution);
double relativeSolutionError(const std::vector<double> &candidate,
                             const std::vector<double> &reference);

} // namespace xfemm_benchmark

#endif
