#include "xfemm_system.h"

#include <iomanip>
#include <iostream>
#include <stdexcept>

int main(int argc, char **argv)
{
    if (argc < 2) {
        std::cerr << "usage: xfemm-system-inspect SYSTEM [SYSTEM ...]\n";
        return 2;
    }
    try {
        std::cout << "path,n,nnz,solve_index,warm_start,tolerance,cpu_threads,"
                     "cpu_iterations,cpu_pcg_residual,true_relative_residual,"
                     "cpu_solution_hash,topology_hash,topology_matches_first\n";
        std::uint64_t firstHash = 0;
        for (int i = 1; i < argc; ++i) {
            const auto system = xfemm_benchmark::readLinearSystem(argv[i]);
            const auto hash = xfemm_benchmark::topologyHash(system);
            if (i == 1) firstHash = hash;
            std::cout << system.sourcePath << ',' << system.dimension() << ','
                      << system.nonzeros() << ',' << system.solveIndex << ','
                      << (system.warmStart() ? "yes" : "no") << ','
                      << std::setprecision(17) << system.tolerance << ','
                      << system.cpuThreads << ',' << system.cpuIterations << ','
                      << system.cpuPcgResidual << ','
                      << xfemm_benchmark::relativeResidual(system,
                                                           system.cpuSolution)
                      << ',' << std::hex
                      << xfemm_benchmark::solutionHash(system.cpuSolution) << ','
                      << hash << std::dec << ','
                      << (hash == firstHash ? "yes" : "no") << '\n';
        }
    } catch (const std::exception &error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
