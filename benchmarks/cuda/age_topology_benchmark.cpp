#include "cudss_solver.h"
#include "xfemm_system.h"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace xb = xfemm_benchmark;

namespace {

void printResult(const char *mode, const xb::BenchmarkResult &result)
{
    std::cout << "system," << mode << ',' << result.system << ','
              << result.dimension << ',' << result.nonzeros << ','
              << std::fixed << std::setprecision(6) << result.setupMs << ','
              << result.uploadMs << ',' << result.analysisMs << ','
              << result.factorizationMs << ',' << result.solveMs << ','
              << result.factorNonzeros << ',' << result.permanentMemoryBytes
              << ',' << result.memoryBytes << ','
              << std::scientific << std::setprecision(9)
              << result.relativeResidual << ',' << result.solutionError << ','
              << result.solutionHash << ',' << result.repeatabilityError << '\n';
}

double sumField(const std::vector<xb::BenchmarkResult> &results,
                double xb::BenchmarkResult::*field)
{
    double sum = 0.;
    for (const auto &result : results) sum += result.*field;
    return sum;
}

} // namespace

int main(int argc, char **argv)
{
    bool completeRevolutionSuperset = true;
    std::string dumpPrefix;
    int firstArgument = 1;
    while (firstArgument < argc && argv[firstArgument][0] == '-') {
        const std::string option = argv[firstArgument++];
        if (option == "--observed-only")
            completeRevolutionSuperset = false;
        else if (option == "--dump-solutions" && firstArgument < argc)
            dumpPrefix = argv[firstArgument++];
        else {
            std::cerr << "unknown or incomplete option: " << option << '\n';
            return 2;
        }
    }
    if (argc - firstArgument < 3) {
        std::cerr << "usage: xfemm-cuda-age-topology-benchmark "
                     "[--observed-only] [--dump-solutions PREFIX] "
                     "REFERENCE SHIFTED SYSTEM [SYSTEM ...]\n";
        return 2;
    }
    try {
        const auto reference = xb::readLinearSystem(argv[firstArgument]);
        const auto shifted = xb::readLinearSystem(argv[firstArgument + 1]);
        std::vector<xb::LinearSystem> systems;
        for (int i = firstArgument + 2; i < argc; ++i)
            systems.push_back(xb::readLinearSystem(argv[i]));

        std::map<std::uint64_t, std::size_t> frequencies;
        for (const auto &system : systems) ++frequencies[xb::topologyHash(system)];
        const auto result =
            xb::runCudssAgeTopologyExperiment(reference, shifted, systems,
                                               completeRevolutionSuperset);

        if (!dumpPrefix.empty()) {
            for (std::size_t i = 0; i < result.universalSolutions.size(); ++i) {
                const std::string path = dumpPrefix + "." + std::to_string(i) + ".bin";
                std::ofstream output(path, std::ios::binary | std::ios::trunc);
                if (!output) throw std::runtime_error("cannot create " + path);
                const std::uint64_t count = result.universalSolutions[i].size();
                output.write(reinterpret_cast<const char *>(&count), sizeof(count));
                output.write(reinterpret_cast<const char *>(
                                 result.universalSolutions[i].data()),
                             static_cast<std::streamsize>(count * sizeof(double)));
                if (!output) throw std::runtime_error("cannot write " + path);
            }
        }

        std::cout << "summary,n,angle_specific_nnz,universal_nnz,nnz_growth_pct,"
                     "changed_edges,ring_a_nodes,ring_b_nodes,changed_components,"
                     "observed_systems,distinct_topologies,union_build_ms,"
                     "universal_lower_nnz,universal_update_bytes,"
                     "specific_cache_resident_bytes,specific_cache_peak_bytes,"
                     "universal_resident_bytes,universal_peak_bytes,"
                     "specific_wall_ms,universal_wall_ms\n";
        std::cout << "summary," << reference.dimension() << ','
                  << result.angleSpecificNonzeros << ',' << result.universalNonzeros
                  << ',' << std::fixed << std::setprecision(6)
                  << (100. * (static_cast<double>(result.universalNonzeros) /
                              result.angleSpecificNonzeros - 1.)) << ','
                  << result.changedUndirectedEdges << ',' << result.firstRingNodes
                  << ',' << result.secondRingNodes << ','
                  << result.changedGraphComponents << ',' << systems.size() << ','
                  << frequencies.size() << ',' << result.universalBuildMs << ','
                  << result.universalLowerNonzeros << ','
                  << result.universalUpdateBytes << ','
                  << result.angleSpecificCacheResidentBytes << ','
                  << result.angleSpecificCachePeakBytes << ','
                  << result.universalResidentBytes << ','
                  << result.universalPeakBytes << ','
                  << result.angleSpecificWallMs << ',' << result.universalWallMs
                  << '\n';

        std::cout << "topology,hash,frequency\n";
        for (const auto &entry : frequencies)
            std::cout << "topology," << std::hex << entry.first << std::dec << ','
                      << entry.second << '\n';

        std::cout << "record,mode,path,n,nnz,setup_ms,upload_ms,analysis_ms,"
                     "factorization_ms,solve_ms,factor_nnz,permanent_memory_bytes,"
                     "peak_memory_bytes,"
                     "relative_residual,solution_error,solution_hash,"
                     "repeatability_error\n";
        for (const auto &item : result.angleSpecific)
            printResult("angle-specific", item);
        for (const auto &item : result.universal)
            printResult("universal-age", item);

        std::cout << "totals,mode,analysis_ms,upload_ms,factorization_ms,solve_ms\n";
        for (const auto *mode : {"angle-specific", "universal-age"}) {
            const auto &items = std::string(mode) == "angle-specific"
                              ? result.angleSpecific : result.universal;
            std::cout << "totals," << mode << ','
                      << sumField(items, &xb::BenchmarkResult::analysisMs) << ','
                      << sumField(items, &xb::BenchmarkResult::uploadMs) << ','
                      << sumField(items, &xb::BenchmarkResult::factorizationMs)
                      << ',' << sumField(items, &xb::BenchmarkResult::solveMs)
                      << '\n';
        }
    } catch (const std::exception &error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
