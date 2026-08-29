#include "PersistentMotorSession.h"

#include "fpproc.h"
#include "femmenums.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>

namespace femm {
namespace {

using Clock = std::chrono::steady_clock;

double elapsedMs(Clock::time_point start)
{
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

bool atLeast(PhysicalResultLevel lhs, PhysicalResultLevel rhs)
{
    return static_cast<int>(lhs) >= static_cast<int>(rhs);
}

double solverElementAreaM2(const FSolver &solver,
                           const femmsolver::CMElement &element)
{
    const auto &n0 = solver.meshnode[static_cast<std::size_t>(element.p[0])];
    const auto &n1 = solver.meshnode[static_cast<std::size_t>(element.p[1])];
    const auto &n2 = solver.meshnode[static_cast<std::size_t>(element.p[2])];
    // FSolver mesh coordinates are centimetres.
    return 0.0001 * std::abs((n1.x - n0.x) * (n2.y - n0.y) -
                             (n2.x - n0.x) * (n1.y - n0.y)) / 2.0;
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
        initializeDirectLinkageExtractor();
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

void PersistentMotorSession::initializeDirectLinkageExtractor()
{
    if (m_directLinkageInitialized)
        return;
    m_directLinkageInitialized = true;

    const auto &problem = m_session.model().problem();
    const auto &solver = m_backend->solvedSolver();
    const std::size_t circuits = problem.circproplist.size();
    if (solver.Frequency != 0 || problem.problemType != PLANAR ||
        circuits != static_cast<std::size_t>(solver.NumCircPropsOrig))
        return;
    for (const auto &property : problem.circproplist) {
        const auto *circuit = dynamic_cast<const CMCircuit *>(property.get());
        if (!circuit || circuit->CircType != 1)
            return;
    }

    std::vector<double> labelArea(solver.labellist.size(), 0.0);
    for (const auto &element : solver.meshele) {
        if (element.lbl >= 0 &&
            element.lbl < static_cast<int>(labelArea.size()))
            labelArea[static_cast<std::size_t>(element.lbl)] +=
                solverElementAreaM2(solver, element);
    }

    const double depth = problem.Depth == -1
        ? 1.0
        : problem.Depth * LengthConvMeters[problem.LengthUnits];
    m_directLinkageWeights.assign(
        circuits, std::vector<double>(solver.meshnode.size(), 0.0));
    std::vector<bool> found(circuits, false);
    for (const auto &element : solver.meshele) {
        if (element.lbl < 0 ||
            element.lbl >= static_cast<int>(solver.labellist.size()))
            continue;
        const std::size_t label = static_cast<std::size_t>(element.lbl);
        const int privateCircuit = solver.labellist[label].InCircuit;
        if (privateCircuit < 0 ||
            privateCircuit >= static_cast<int>(solver.circproplist.size()))
            continue;
        const int originalCircuit = solver.circproplist[
            static_cast<std::size_t>(privateCircuit)].OrigCirc;
        if (originalCircuit < 0 ||
            originalCircuit >= static_cast<int>(circuits) ||
            !(labelArea[label] > 0))
            continue;
        found[static_cast<std::size_t>(originalCircuit)] = true;
        const double coefficient =
            static_cast<double>(solver.labellist[label].Turns) * depth *
            solverElementAreaM2(solver, element) / (3.0 * labelArea[label]);
        auto &weights = m_directLinkageWeights[
            static_cast<std::size_t>(originalCircuit)];
        for (int node : element.p)
            weights[static_cast<std::size_t>(node)] += coefficient;
    }
    if (!std::all_of(found.begin(), found.end(), [](bool value) { return value; })) {
        m_directLinkageWeights.clear();
        return;
    }
    m_directLinkageAvailable = true;
}

void PersistentMotorSession::refreshPhysicalView(bool airGapOnly)
{
    if (!m_physicalViewInitialized) {
        initializePhysicalView();
        return;
    }
    const bool updated = airGapOnly
        ? m_postProcessor->UpdateAirGapSolution(m_backend->solvedSolver(),
                                                m_backend->solvedSystem())
        : m_postProcessor->UpdateSolution(m_backend->solvedSolver(),
                                          m_backend->solvedSystem());
    if (!updated)
        throw std::runtime_error("could not refresh persistent magnetic physical view");
}

TrialSolution PersistentMotorSession::evaluateTrial(PhysicalResultLevel level)
{
    const auto totalStarted = Clock::now();
    TrialSolution trial = m_session.solve();
    m_lastEvaluatedTrialId = trial.id;
    if (m_initialized)
        trial.diagnostics.sessionInitializationMs =
            m_initializationDiagnostics.sessionInitializationMs;
    initializeDirectLinkageExtractor();

    const bool needsFullView =
        level == PhysicalResultLevel::FullDiagnostics ||
        !m_directLinkageAvailable;
    if (needsFullView) {
        const auto updateStarted = Clock::now();
        refreshPhysicalView(false);
        trial.diagnostics.serializationPostprocessorMs += elapsedMs(updateStarted);
    }
    addLinkageResults(trial, level == PhysicalResultLevel::FullDiagnostics);
    if (atLeast(level, PhysicalResultLevel::AcceptedState)) {
        if (!needsFullView) {
            const auto updateStarted = Clock::now();
            refreshPhysicalView(true);
            trial.diagnostics.serializationPostprocessorMs += elapsedMs(updateStarted);
        }
        addAcceptedStateResults(
            trial, level == PhysicalResultLevel::FullDiagnostics);
    }
    if (level == PhysicalResultLevel::FullDiagnostics)
        addFullDiagnosticResults(trial);
    trial.real->physicalResultLevel = level;
    trial.diagnostics.totalEvaluateMs = elapsedMs(totalStarted);
    return trial;
}

void PersistentMotorSession::completeTrial(TrialSolution &trial,
                                           PhysicalResultLevel level)
{
    if (!trial.real)
        throw std::logic_error("persistent motor evaluator requires a real solution");
    if (trial.id != m_lastEvaluatedTrialId)
        throw std::logic_error("only the latest live magnetic trial can be completed");
    const PhysicalResultLevel current = trial.real->physicalResultLevel;
    if (atLeast(current, level))
        return;

    const auto completionStarted = Clock::now();
    const auto updateStarted = Clock::now();
    refreshPhysicalView(level == PhysicalResultLevel::AcceptedState);
    trial.diagnostics.serializationPostprocessorMs += elapsedMs(updateStarted);
    if (level == PhysicalResultLevel::FullDiagnostics) {
        addLinkageResults(trial, true);
        addAcceptedStateResults(trial, true);
        addFullDiagnosticResults(trial);
    } else {
        addAcceptedStateResults(trial, false);
    }
    trial.real->physicalResultLevel = level;
    trial.diagnostics.totalEvaluateMs += elapsedMs(completionStarted);
}

TrialSolution PersistentMotorSession::evaluateTrial(
    double thetaMechanical, const std::array<double, 3> &currents,
    PhysicalResultLevel level)
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
    return evaluateTrial(level);
}

void PersistentMotorSession::addLinkageResults(
    TrialSolution &trial, bool validateWithPostProcessor)
{
    if (!trial.real)
        throw std::logic_error("persistent motor evaluator requires a real solution");
    if (m_directLinkageAvailable &&
        m_directLinkageWeights.size() != trial.real->circuits.size())
        throw std::logic_error("direct linkage circuit count changed");

    const auto fluxStarted = Clock::now();
    for (std::size_t i = 0; i < trial.real->circuits.size(); ++i) {
        double linkage = 0;
        if (m_directLinkageAvailable) {
            const auto &weights = m_directLinkageWeights[i];
            const auto &solution = m_backend->solvedSystem().rhs();
            if (weights.size() != solution.size())
                throw std::logic_error("direct linkage solution dimension changed");
            for (std::size_t node = 0; node < weights.size(); ++node)
                linkage += weights[node] * solution[node];
            trial.real->circuits[i].strandedFluxLinkage = linkage;
        } else {
            const CComplex conventional = m_postProcessor->GetFluxLinkage(
                static_cast<int>(i));
            const auto stranded = strandedSeriesFluxLinkage(
                *m_postProcessor, static_cast<int>(i));
            linkage = stranded.value_or(conventional).re;
            trial.real->circuits[i].conventionalFluxLinkage = conventional.re;
            if (stranded)
                trial.real->circuits[i].strandedFluxLinkage = stranded->re;
        }
        if (validateWithPostProcessor) {
            const CComplex conventional = m_postProcessor->GetFluxLinkage(
                static_cast<int>(i));
            trial.real->circuits[i].conventionalFluxLinkage = conventional.re;
            const auto reference = strandedSeriesFluxLinkage(
                *m_postProcessor, static_cast<int>(i));
            if (reference) {
                const double tolerance = 1e-13 + 1e-10 *
                    std::max(std::abs(linkage), std::abs(reference->re));
                if (std::abs(linkage - reference->re) > tolerance)
                    throw std::runtime_error(
                        "direct stranded-series linkage disagrees with FPProc");
            }
        }
        trial.real->circuits[i].fluxLinkage = linkage;
        trial.circuits[i].fluxLinkage = CComplex(linkage, 0);
    }
    trial.diagnostics.fluxLinkageMs += elapsedMs(fluxStarted);
}

void PersistentMotorSession::addAcceptedStateResults(
    TrialSolution &trial, bool packageHarmonics)
{
    if (!trial.real)
        throw std::logic_error("persistent motor evaluator requires a real solution");
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
        if (packageHarmonics) {
            const auto harmonicStarted = Clock::now();
            result.harmonics.reserve(static_cast<std::size_t>(gap.nn));
            for (int i = 0; i < gap.nn; ++i)
                result.harmonics.push_back({gap.nh[i], gap.brc[i], gap.brs[i],
                                            gap.btc[i], gap.bts[i]});
            trial.diagnostics.airGapHarmonicPackagingMs +=
                elapsedMs(harmonicStarted);
        }
        trial.real->torque += result.torque;
        trial.real->airGaps.push_back(std::move(result));
    }
}

void PersistentMotorSession::addFullDiagnosticResults(TrialSolution &trial)
{
    if (!trial.real)
        throw std::logic_error("persistent motor evaluator requires a real solution");
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
    trial.diagnostics.energyCoenergyMs += elapsedMs(energyStarted);
}

std::shared_ptr<const AcceptedState>
PersistentMotorSession::commitTrial(const TrialSolution &trial)
{
    return m_session.acceptSolution(trial);
}

void PersistentMotorSession::rollbackToCommitted()
{
    m_session.rollbackToCommitted();
    // The backend owns the authoritative committed/trial magnetic state. The
    // physical view is refreshed explicitly by the next accepted/full result.
}

} // namespace femm
