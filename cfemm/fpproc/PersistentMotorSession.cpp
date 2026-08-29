#include "PersistentMotorSession.h"

#include "fpproc.h"

#include <chrono>
#include <stdexcept>

namespace femm {
namespace {

using Clock = std::chrono::steady_clock;

double elapsedMs(Clock::time_point start)
{
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

std::optional<CComplex> strandedSeriesFluxLinkage(
    const FPProc &postProcessor, int circuit)
{
    // GetFluxLinkage's A.J/conj(I) path is algebraically equivalent for a
    // nonzero series-circuit current, but it changes algorithms at I == 0.
    // Transient finite differences must remain continuous through current
    // zero.  For DC planar stranded windings, directly summing turns times
    // average A is valid at every current and is already the established FEMM
    // zero-current calculation.
    if (postProcessor.Frequency == 0 &&
        postProcessor.problemType == femm::PLANAR &&
        postProcessor.circproplist.at(static_cast<std::size_t>(circuit)).CircType == 1) {
        CComplex linkage;
        for (std::size_t label = 0; label < postProcessor.blocklist.size(); ++label)
            if (postProcessor.blocklist[label].InCircuit == circuit)
                linkage += postProcessor.GetStrandedLinkage(
                    static_cast<int>(label));
        return linkage;
    }
    return std::nullopt;
}

} // namespace

PersistentMotorSession::PersistentMotorSession(ModelDefinition model)
    : m_backend(std::make_shared<FSolverAnalysisBackend>()),
      m_session(std::move(model), m_backend),
      m_postProcessor(std::make_unique<FPProc>())
{}

PersistentMotorSession::PersistentMotorSession(ModelDefinition model,
                                               CudssSessionOptions options)
    : m_backend(std::make_shared<FSolverAnalysisBackend>(options)),
      m_session(std::move(model), m_backend),
      m_postProcessor(std::make_unique<FPProc>())
{}

PersistentMotorSession::~PersistentMotorSession() = default;

void PersistentMotorSession::initialize()
{
    if (m_initialized)
        return;
    const auto started = Clock::now();
    m_session.initialize();
    m_initializationDiagnostics = m_backend->initializationDiagnostics();
    if (m_backend->solveCount() != 0) {
        const auto physicalViewStarted = Clock::now();
        initializePhysicalView();
        m_initializationDiagnostics.serializationPostprocessorMs +=
            elapsedMs(physicalViewStarted);
    }
    m_initializationDiagnostics.sessionInitializationMs = elapsedMs(started);
    m_initializationDiagnostics.totalEvaluateMs =
        m_initializationDiagnostics.sessionInitializationMs;
    m_initialized = true;
}

void PersistentMotorSession::initializePhysicalView()
{
    if (m_physicalViewInitialized)
        return;
    if (!m_postProcessor->OpenDocument(m_session.model().problem(),
                                       m_backend->solvedSolver(),
                                       m_backend->solvedSystem()))
        throw std::runtime_error("could not initialize persistent magnetic physical view");
    m_physicalViewInitialized = true;
}

TrialSolution PersistentMotorSession::evaluateTrial()
{
    const auto totalStarted = Clock::now();
    TrialSolution trial = m_session.solve();
    if (m_initialized)
        trial.diagnostics.sessionInitializationMs =
            m_initializationDiagnostics.sessionInitializationMs;
    const auto updateStarted = Clock::now();
    if (!m_physicalViewInitialized)
        initializePhysicalView();
    else if (!m_postProcessor->UpdateSolution(m_backend->solvedSolver(),
                                               m_backend->solvedSystem()))
        throw std::runtime_error("could not refresh persistent magnetic physical view");
    trial.diagnostics.serializationPostprocessorMs += elapsedMs(updateStarted);
    addPhysicalResults(trial);
    trial.diagnostics.totalEvaluateMs = elapsedMs(totalStarted);
    return trial;
}

TrialSolution PersistentMotorSession::evaluateTrial(
    double thetaMechanical, const std::array<double, 3> &currents)
{
    if (m_session.model().problem().circproplist.size() != currents.size())
        throw std::invalid_argument(
            "three-phase evaluateTrial requires exactly three model circuits");
    m_session.updateSolveParameters([&](SolveParameters &parameters) {
        for (std::size_t i = 0; i < currents.size(); ++i)
            parameters.circuitConstraints[CircuitId{i}] = {
                CircuitConstraintKind::PrescribedCurrent, CComplex(currents[i], 0)};
        for (auto &gap : parameters.airGapPositions)
            gap.second.innerAngle = thetaMechanical;
    });
    return evaluateTrial();
}

void PersistentMotorSession::addPhysicalResults(TrialSolution &trial)
{
    if (!trial.real)
        throw std::logic_error("persistent motor evaluator requires a real solution");

    const auto fluxStarted = Clock::now();
    for (std::size_t i = 0; i < trial.real->circuits.size(); ++i) {
        const CComplex conventional = m_postProcessor->GetFluxLinkage(
            static_cast<int>(i));
        const auto stranded = strandedSeriesFluxLinkage(
            *m_postProcessor, static_cast<int>(i));
        const CComplex linkage = stranded.value_or(conventional);
        trial.real->circuits[i].fluxLinkage = linkage.re;
        trial.real->circuits[i].conventionalFluxLinkage = conventional.re;
        if (stranded)
            trial.real->circuits[i].strandedFluxLinkage = stranded->re;
        trial.circuits[i].fluxLinkage = linkage;
    }
    trial.diagnostics.fluxLinkageMs = elapsedMs(fluxStarted);

    trial.real->airGaps.clear();
    trial.real->torque = 0;
    trial.real->airGaps.reserve(m_postProcessor->agelist.size());
    for (const auto &gap : m_postProcessor->agelist) {
        RealAirGapResult result;
        result.name = gap.BdryName;
        result.centerVectorPotential = gap.aco;
        const auto torqueStarted = Clock::now();
        if (m_postProcessor->gapDCTorqueIntegral(gap.BdryName, result.torque) !=
            FPProcError::NoError)
            throw std::runtime_error("could not calculate AGE torque for " + gap.BdryName);
        trial.diagnostics.torqueMs += elapsedMs(torqueStarted);
        const auto harmonicStarted = Clock::now();
        result.harmonics.reserve(static_cast<std::size_t>(gap.nn));
        for (int i = 0; i < gap.nn; ++i)
            result.harmonics.push_back({gap.nh[i], gap.brc[i], gap.brs[i],
                                        gap.btc[i], gap.bts[i]});
        trial.diagnostics.airGapHarmonicPackagingMs +=
            elapsedMs(harmonicStarted);
        trial.real->torque += result.torque;
        trial.real->airGaps.push_back(std::move(result));
    }

    const auto energyStarted = Clock::now();
    std::vector<bool> selected;
    selected.reserve(m_postProcessor->blocklist.size());
    for (auto &label : m_postProcessor->blocklist) {
        selected.push_back(label.IsSelected);
        label.IsSelected = true;
    }
    trial.real->magneticFieldEnergyJ = m_postProcessor->BlockIntegral(2).re;
    trial.real->magneticFieldCoenergyJ = m_postProcessor->BlockIntegral(17).re;
    for (const auto &gap : m_postProcessor->agelist) {
        CComplex gapEnergy;
        if (m_postProcessor->gapTimeAvgStoredEnergyIntegral(
                gap.BdryName, gapEnergy) != FPProcError::NoError)
            throw std::runtime_error(
                "could not calculate AGE stored energy for " + gap.BdryName);
        // The AGE is a linear air region, so its energy and coenergy are equal.
        trial.real->magneticFieldEnergyJ += gapEnergy.re;
        trial.real->magneticFieldCoenergyJ += gapEnergy.re;
    }
    for (std::size_t i = 0; i < selected.size(); ++i)
        m_postProcessor->blocklist[i].IsSelected = selected[i];
    trial.diagnostics.energyCoenergyMs = elapsedMs(energyStarted);
}

std::shared_ptr<const AcceptedState>
PersistentMotorSession::commitTrial(const TrialSolution &trial)
{
    return m_session.acceptSolution(trial);
}

void PersistentMotorSession::rollbackToCommitted()
{
    m_session.rollbackToCommitted();
    if (m_physicalViewInitialized &&
        !m_postProcessor->UpdateSolution(m_backend->solvedSolver(),
                                         m_backend->solvedSystem()))
        throw std::runtime_error("could not restore committed physical state");
}

} // namespace femm
