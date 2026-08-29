#ifndef XFEMM_PERSISTENT_MOTOR_SESSION_H
#define XFEMM_PERSISTENT_MOTOR_SESSION_H

#include "AnalysisSession.h"
#include "FSolverAnalysisBackend.h"

#include <array>
#include <cstdint>
#include <memory>
#include <vector>

class FPProc;

namespace femm {

/**
 * Persistent, in-memory magnetic evaluator.
 *
 * The model, mesh, ordering, FSolver, linear backend, cuDSS buckets and one
 * physical-calculation view all outlive individual evaluateTrial() calls.
 * The default constructor is the unchanged legacy CPU reference. Supplying
 * CudssSessionOptions opts into the explicit overload.
 */
class PersistentMotorSession {
public:
    explicit PersistentMotorSession(ModelDefinition model);
    PersistentMotorSession(ModelDefinition model, CudssSessionOptions options);
    ~PersistentMotorSession();

    PersistentMotorSession(const PersistentMotorSession &) = delete;
    PersistentMotorSession &operator=(const PersistentMotorSession &) = delete;

    AnalysisSession &analysis() { return m_session; }
    const AnalysisSession &analysis() const { return m_session; }
    std::size_t topologyImportCount() const { return m_backend->topologyImportCount(); }
    std::size_t orderingCount() const { return m_backend->orderingCount(); }
    std::size_t meshFileReadCount() const { return m_backend->meshFileReadCount(); }
    std::size_t meshFileWriteCount() const { return m_backend->meshFileWriteCount(); }
    std::size_t residentBucketCount() const { return m_backend->residentBucketCount(); }
    std::size_t bucketEvictionCount() const { return m_backend->bucketEvictionCount(); }
    std::size_t bucketDefinitionCount() const { return m_backend->bucketDefinitionCount(); }

    /**
     * Pay the initial cuDSS bucket construction/symbolic/factorization and
     * initialize the reusable physical view before entering the hot loop.
     */
    void initialize();
    const EvaluationDiagnostics &initializationDiagnostics() const
    { return m_initializationDiagnostics; }

    /**
     * Evaluate the currently configured angle/circuit state as a trial.
     * FullDiagnostics is the compatibility default; nonlinear electrical
     * solvers should request ResidualOnly and promote only converged states.
     */
    TrialSolution evaluateTrial(
        PhysicalResultLevel level = PhysicalResultLevel::FullDiagnostics);

    /** Add accepted-state or full diagnostics to the latest live trial. */
    void completeTrial(TrialSolution &trial, PhysicalResultLevel level);

    /**
     * Three-phase convenience interface moving toward
     * evaluate(theta_mech,currents)->lambda[3],torque,diagnostics.
     * The model must contain exactly three circuits.
     */
    TrialSolution evaluateTrial(double thetaMechanical,
                                const std::array<double, 3> &currents,
                                PhysicalResultLevel level =
                                    PhysicalResultLevel::FullDiagnostics);

    std::shared_ptr<const AcceptedState> commitTrial(const TrialSolution &trial);
    void rollbackToCommitted();

private:
    void initializePhysicalView();
    void initializeDirectLinkageExtractor();
    void refreshPhysicalView(bool airGapOnly);
    void addLinkageResults(TrialSolution &trial, bool validateWithPostProcessor);
    void addAcceptedStateResults(TrialSolution &trial, bool packageHarmonics);
    void addFullDiagnosticResults(TrialSolution &trial);

    std::shared_ptr<FSolverAnalysisBackend> m_backend;
    AnalysisSession m_session;
    std::unique_ptr<FPProc> m_postProcessor;
    bool m_physicalViewInitialized = false;
    bool m_directLinkageInitialized = false;
    bool m_directLinkageAvailable = false;
    bool m_initialized = false;
    std::uint64_t m_lastEvaluatedTrialId = 0;
    std::vector<std::vector<double>> m_directLinkageWeights;
    EvaluationDiagnostics m_initializationDiagnostics;
};

} // namespace femm

#endif
