#ifndef XFEMM_FSOLVERANALYSISBACKEND_H
#define XFEMM_FSOLVERANALYSISBACKEND_H

#include "AnalysisSession.h"
#include "fsolver.h"

#include <cstddef>
#include <map>
#include <memory>
#include <string>
#include <vector>

#ifdef XFEMM_USE_CUDSS
#include "linsolve/CudssLinearSystemBackend.h"
#endif

namespace femm {

/** AnalysisSession adapter for the legacy magnetic FSolver engine.
 *
 * The imported mesh and its Cuthill--McKee numbering are owned by this
 * backend and replaced only when meshTopologyIdentity changes. The default
 * path creates the unchanged Legacy system per solve; the opt-in cuDSS path
 * retains bucket analysis, factorization workspace, and the final tangent.
 */
struct CudssSessionOptions {
    /** Conservative width selected by the measured continuous AGE experiment. */
    double bucketWidthDegrees = 1.2;
    bool deterministic = false;
    /** If true, a failed CUDA evaluation is restarted explicitly on Legacy. */
    bool explicitLegacyFallback = false;
    /** Bounded resident cuDSS bucket cache; two retains committed plus trial. */
    std::size_t bucketCacheCapacity = 2;
};

class FSolverAnalysisBackend final : public AnalysisSolverBackend {
public:
    /** Unchanged legacy CPU reference path. */
    FSolverAnalysisBackend();
    /** Opt-in persistent native-real FP64 cuDSS path. */
    explicit FSolverAnalysisBackend(CudssSessionOptions options);
    ~FSolverAnalysisBackend() override;

    void synchronize(const ModelDefinition &, const SolveParameters &,
                     const PreparedAnalysis &, std::shared_ptr<const mesh::SolverMesh>,
                     std::uint64_t meshTopologyIdentity, Dirty rebuilt) override;
    TrialSolution solve(const ModelDefinition &, const SolveParameters &,
                        const PreparedAnalysis &) override;
    void initialize(const ModelDefinition &, const SolveParameters &,
                    const PreparedAnalysis &) override;
    void commitTrial(const TrialSolution &) override;
    void rollbackToCommitted() override;

    /** Write the most recently solved field in the legacy .ans format. */
    void writeSolution(const std::string &ansPath);

    /** Native solved state used to construct an in-memory post-processor view. */
    const FSolver &solvedSolver() const;
    const femm::LinearSystemBackend<double> &solvedSystem() const;

    std::size_t topologyImportCount() const { return m_topologyImports; }
    std::size_t orderingCount() const { return m_orderings; }
    std::size_t couplingRegenerationCount() const { return m_couplingRegenerations; }
    std::size_t operatorAssemblyCount() const { return m_operatorAssemblies; }
    std::size_t rightHandSideAssemblyCount() const { return m_rightHandSideAssemblies; }
    std::size_t solveCount() const { return m_solves; }
    std::size_t meshFileReadCount() const { return m_meshFileReads; }
    std::size_t meshFileWriteCount() const { return m_meshFileWrites; }
    bool usesPersistentCudss() const { return m_useCudss; }
    const EvaluationDiagnostics &initializationDiagnostics() const
    { return m_initializationDiagnostics; }
    std::size_t residentBucketCount() const;
    std::size_t bucketEvictionCount() const;
    std::size_t bucketDefinitionCount() const
    {
#ifdef XFEMM_USE_CUDSS
        return m_cudssBuckets.size();
#else
        return 0;
#endif
    }

private:
    void configure(const ModelDefinition &, const SolveParameters &,
                   const PreparedAnalysis &);
    void positionAirGaps(const PreparedAnalysis &);
    TrialSolution solveConfigured(const ModelDefinition &, const SolveParameters &,
                                  const PreparedAnalysis &, bool initialization);
#ifdef XFEMM_USE_CUDSS
    std::string cudssBucketIdentity(const PreparedAnalysis &) const;
    CudssBucketDefinition buildCudssBucket(const PreparedAnalysis &) const;
#endif
    std::unique_ptr<FSolver> m_solver;
    std::unique_ptr<femm::LinearSystemBackend<double>> m_lastSystem;
    std::shared_ptr<const mesh::SolverMesh> m_mesh;
    std::uint64_t m_topologyIdentity = 0;
    bool m_useCudss = false;
    CudssSessionOptions m_cudssOptions;
#ifdef XFEMM_USE_CUDSS
    CudssBucketDefinition m_currentBucket;
    std::string m_committedBucketIdentity;
    std::map<std::string, CudssBucketDefinition> m_cudssBuckets;
#endif
    std::vector<double> m_committedSolution;
    std::vector<double> m_trialSolution;
    bool m_haveCommittedSolution = false;
    bool m_haveTrialSolution = false;
    double m_lastAirGapUpdateMs = 0;
    double m_pendingBucketConstructionMs = 0;
    double m_initializationMs = 0;
    EvaluationDiagnostics m_initializationDiagnostics;
    bool m_initialized = false;
    std::size_t m_topologyImports = 0;
    std::size_t m_orderings = 0;
    std::size_t m_couplingRegenerations = 0;
    std::size_t m_operatorAssemblies = 0;
    std::size_t m_rightHandSideAssemblies = 0;
    std::size_t m_solves = 0;
    // The session path is in-memory. These counters make accidental legacy
    // mesh-file I/O observable to callers and regression tests.
    std::size_t m_meshFileReads = 0;
    std::size_t m_meshFileWrites = 0;
};

} // namespace femm

#endif
