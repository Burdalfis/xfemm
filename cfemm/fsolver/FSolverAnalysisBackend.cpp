#include "FSolverAnalysisBackend.h"

#include "CBoundaryProp.h"
#include "CBlockLabel.h"
#include "CCircuit.h"
#include "CMaterialProp.h"
#include "CPointProp.h"
#include "linsolve/backend_factory.h"
#include "femmconstants.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace femm {
namespace {

template<class Target, class Source>
Target magneticCopy(const std::unique_ptr<Source> &source, const char *description)
{
    const auto *typed = dynamic_cast<const Target *>(source.get());
    if (!typed)
        throw std::invalid_argument(std::string("magnetic model contains a non-magnetic ") + description);
    return *typed;
}

double canonicalDegrees(double angle)
{
    double result = std::fmod(angle, 360.0);
    if (result < 0) result += 360.0;
    return result;
}

void positionGap(femmsolver::CAirGapElement &gap, double innerAngle,
                 double outerAngle)
{
    if (gap.innerRingTopology.empty() || gap.outerRingTopology.empty()) {
        gap.InnerAngle = innerAngle;
        gap.OuterAngle = outerAngle;
        return;
    }
    const double step = gap.totalArcLength / gap.totalArcElements;
    auto positioned = [step](const std::vector<CQuadPoint> &topology, double angle) {
        auto ring = topology;
        for (auto &point : ring) {
            point.w0 = canonicalDegrees(point.w0 * step + angle) / step;
        }
        std::stable_sort(ring.begin(), ring.end(), [](const CQuadPoint &a,
                                                      const CQuadPoint &b) {
            return a.w0 < b.w0;
        });
        return ring;
    };
    const auto inner = positioned(gap.innerRingTopology, innerAngle);
    const auto outer = positioned(gap.outerRingTopology, outerAngle);
    const int fullCount = static_cast<int>(inner.size());
    if (fullCount == 0 || outer.size() != inner.size())
        throw std::runtime_error("invalid persistent AGE ring topology");
    gap.InnerShift = inner.front().w0;
    gap.OuterShift = outer.front().w0;
    gap.quadNode.clear();
    gap.quadNode.reserve(gap.totalArcElements + 1);
    for (int i = 0; i <= gap.totalArcElements; ++i) {
        const int p1 = i == fullCount ? 0 : i;
        const int p0 = p1 == 0 ? fullCount - 1 : p1 - 1;
        CQuadPoint q;
        q.n0=inner[p0].n0; q.n1=inner[p1].n0;
        q.n2=outer[p0].n0; q.n3=outer[p1].n0;
        q.w0=inner[p0].w1; q.w1=inner[p1].w1;
        q.w2=outer[p0].w1; q.w3=outer[p1].w1;
        gap.quadNode.push_back(q);
    }
    gap.InnerAngle = canonicalDegrees(innerAngle);
    gap.OuterAngle = canonicalDegrees(outerAngle);
}

#ifdef XFEMM_USE_CUDSS
void appendAgeStructure(
    const femmsolver::CAirGapElement &gap,
    std::vector<std::pair<std::int32_t, std::int32_t>> &entries)
{
    int nodes[10];
    for (int k = 0; k < gap.totalArcElements; ++k) {
        nodes[0] = gap.quadNode[k == 0 ? gap.totalArcElements - 1 : k - 1].n0;
        nodes[1] = gap.quadNode[k].n0;
        nodes[2] = gap.quadNode[k].n1;
        nodes[3] = gap.quadNode[k + 1].n1;
        nodes[4] = gap.quadNode[(k + 2) > gap.totalArcElements ? 1 : k + 2].n1;
        nodes[5] = gap.quadNode[k == 0 ? gap.totalArcElements - 1 : k - 1].n2;
        nodes[6] = gap.quadNode[k].n2;
        nodes[7] = gap.quadNode[k].n3;
        nodes[8] = gap.quadNode[k + 1].n3;
        nodes[9] = gap.quadNode[(k + 2) > gap.totalArcElements ? 1 : k + 2].n3;
        for (int i = 0; i < 10; ++i)
            for (int j = i; j < 10; ++j) {
                std::int32_t row = nodes[i], column = nodes[j];
                if (column < row) std::swap(row, column);
                entries.emplace_back(row, column);
            }
    }
}

void closeStructureUnderPeriodicConstraints(
    const FSolver &solver,
    std::vector<std::pair<std::int32_t, std::int32_t>> &entries)
{
    std::vector<int> parent(static_cast<std::size_t>(solver.NumNodes));
    for (int i = 0; i < solver.NumNodes; ++i) parent[static_cast<std::size_t>(i)] = i;
    const auto root = [&parent](int node) {
        int result = node;
        while (parent[static_cast<std::size_t>(result)] != result)
            result = parent[static_cast<std::size_t>(result)];
        return result;
    };
    for (const auto &constraint : solver.pbclist) {
        const int a = root(constraint.x);
        const int b = root(constraint.y);
        if (a != b) parent[static_cast<std::size_t>(b)] = a;
    }

    std::vector<std::vector<int>> equivalent(static_cast<std::size_t>(solver.NumNodes));
    for (int node = 0; node < solver.NumNodes; ++node)
        equivalent[static_cast<std::size_t>(root(node))].push_back(node);

    const auto source = entries;
    for (const auto &entry : source) {
        const auto &rows = equivalent[static_cast<std::size_t>(root(entry.first))];
        const auto &columns = equivalent[static_cast<std::size_t>(root(entry.second))];
        for (int row : rows)
            for (int column : columns) {
                std::int32_t first = row;
                std::int32_t second = column;
                if (second < first) std::swap(first, second);
                entries.emplace_back(first, second);
            }
    }
}
#endif

} // namespace

FSolverAnalysisBackend::FSolverAnalysisBackend() : m_solver(new FSolver) {}

FSolverAnalysisBackend::FSolverAnalysisBackend(CudssSessionOptions options)
    : m_solver(new FSolver), m_useCudss(true), m_cudssOptions(options)
{
    if (!std::isfinite(options.bucketWidthDegrees) ||
        options.bucketWidthDegrees <= 0 || options.bucketWidthDegrees > 360)
        throw std::invalid_argument("cuDSS AGE bucket width must be in (0, 360] degrees");
#ifndef XFEMM_USE_CUDSS
    throw std::runtime_error(
        "persistent cuDSS support was not built; configure with -DXFEMM_USE_CUDSS=ON");
#endif
}

FSolverAnalysisBackend::~FSolverAnalysisBackend() = default;

void FSolverAnalysisBackend::configure(const ModelDefinition &model,
                                       const SolveParameters &parameters,
                                       const PreparedAnalysis &prepared)
{
    const auto &problem = model.problem();
    m_solver->FileFormat = problem.FileFormat;
    m_solver->Frequency = parameters.frequency;
    m_solver->Precision = problem.Precision;
    m_solver->MinAngle = problem.MinAngle;
    m_solver->Depth = problem.Depth;
    m_solver->LengthUnits = problem.LengthUnits;
    m_solver->Coords = problem.Coords;
    m_solver->ProblemType = problem.problemType;
    m_solver->extZo = problem.extZo; m_solver->extRo = problem.extRo; m_solver->extRi = problem.extRi;
    m_solver->ACSolver = problem.ACSolver;
    m_solver->PrevType = problem.PrevType;
    m_solver->previousSolutionFile = problem.previousSolutionFile;
    m_solver->Relax = 1.0;

    m_solver->nodeproplist.clear();
    for (const auto &item : problem.nodeproplist)
        m_solver->nodeproplist.push_back(magneticCopy<CMPointProp>(item, "point property"));
    m_solver->lineproplist.clear();
    for (const auto &item : problem.lineproplist)
        m_solver->lineproplist.push_back(magneticCopy<CMBoundaryProp>(item, "boundary property"));
    m_solver->blockproplist.clear();
    for (const auto &item : prepared.materials) {
        CMSolverMaterialProp solverMaterial;
        static_cast<CMMaterialProp &>(solverMaterial) = item;
        m_solver->blockproplist.push_back(std::move(solverMaterial));
    }
    m_solver->labellist.clear();
    for (const auto &item : problem.labellist)
        m_solver->labellist.push_back(magneticCopy<CMBlockLabel>(item, "block label"));
    m_solver->circproplist.clear();
    for (std::size_t i = 0; i < problem.circproplist.size(); ++i) {
        auto circuit = magneticCopy<CMCircuit>(problem.circproplist[i], "circuit");
        const auto constraint = parameters.circuitConstraints.at(CircuitId{i});
        circuit.Amps = constraint.value;
        circuit.Case = constraint.kind == CircuitConstraintKind::PrescribedCurrent ? 0 : 2;
        circuit.dVolts = constraint.kind == CircuitConstraintKind::PrescribedVoltage
                       ? constraint.value : CComplex();
        m_solver->circproplist.push_back(circuit);
    }

    // Match FSolver::LoadProblemFile's series-circuit preprocessing.  Each
    // physical stranded block gets a private flat-current-density circuit so
    // its signed Turns multiplier is applied.  The persistent in-memory path
    // does not call LoadProblemFile and must perform this once-per-configure
    // transformation explicitly.
    const std::size_t originalCircuitCount = m_solver->circproplist.size();
    for (std::size_t i = 0; i < originalCircuitCount; ++i)
        m_solver->circproplist[i].OrigCirc = -1;
    for (auto &label : m_solver->labellist) {
        if (label.InCircuit < 0)
            continue;
        const std::size_t original = static_cast<std::size_t>(label.InCircuit);
        if (original >= originalCircuitCount)
            throw std::runtime_error("block label references an invalid circuit");
        if (m_solver->circproplist[original].CircType != 1)
            continue;
        CMCircuit blockCircuit = m_solver->circproplist[original];
        blockCircuit.OrigCirc = static_cast<int>(original);
        blockCircuit.Amps *= static_cast<double>(label.Turns);
        label.InCircuit = static_cast<int>(m_solver->circproplist.size());
        m_solver->circproplist.push_back(std::move(blockCircuit));
    }
    for (auto &circuit : m_solver->circproplist)
        if (circuit.CircType == 1)
            circuit.CircType = 0;
    m_solver->NumPointProps = static_cast<int>(m_solver->nodeproplist.size());
    m_solver->NumLineProps = static_cast<int>(m_solver->lineproplist.size());
    m_solver->NumBlockProps = static_cast<int>(m_solver->blockproplist.size());
    m_solver->NumBlockLabels = static_cast<int>(m_solver->labellist.size());
    m_solver->NumCircPropsOrig = static_cast<int>(originalCircuitCount);
    m_solver->NumCircProps = static_cast<int>(m_solver->circproplist.size());
}

void FSolverAnalysisBackend::positionAirGaps(const PreparedAnalysis &prepared)
{
    const auto started = std::chrono::steady_clock::now();
    for (const auto &entry : prepared.airGapPositions)
        for (auto &gap : m_solver->agelist) {
            if (gap.BdryName != m_solver->lineproplist[entry.first.value].BdryName) continue;
            positionGap(gap, entry.second.innerAngle, entry.second.outerAngle);
            ++m_couplingRegenerations;
        }
    m_lastAirGapUpdateMs = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - started).count();
}

#ifdef XFEMM_USE_CUDSS
std::string
FSolverAnalysisBackend::cudssBucketIdentity(const PreparedAnalysis &prepared) const
{
    std::ostringstream identity;
    identity.precision(17);
    identity << "age-width=" << m_cudssOptions.bucketWidthDegrees;
    for (const auto &entry : prepared.airGapPositions) {
        const double relative = canonicalDegrees(entry.second.innerAngle -
                                                 entry.second.outerAngle);
        const int index = static_cast<int>(
            std::floor(relative / m_cudssOptions.bucketWidthDegrees));
        identity << '|' << entry.first.value << ':' << index;
    }
    return identity.str();
}

CudssBucketDefinition
FSolverAnalysisBackend::buildCudssBucket(const PreparedAnalysis &prepared) const
{
    CudssBucketDefinition definition;
    definition.identity = cudssBucketIdentity(prepared);
    for (const auto &entry : prepared.airGapPositions) {
        const std::string &name = m_solver->lineproplist[entry.first.value].BdryName;
        auto found = std::find_if(m_solver->agelist.begin(), m_solver->agelist.end(),
            [&](const femmsolver::CAirGapElement &gap) { return gap.BdryName == name; });
        if (found == m_solver->agelist.end()) continue;

        const double relative = canonicalDegrees(entry.second.innerAngle -
                                                 entry.second.outerAngle);
        const int index = static_cast<int>(
            std::floor(relative / m_cudssOptions.bucketWidthDegrees));
        const double begin = index * m_cudssOptions.bucketWidthDegrees;
        const double end = std::min(360.0, begin + m_cudssOptions.bucketWidthDegrees);

        // AGE topology events occur at ring-point crossings. Sample both
        // sides of every possible event in this narrow interval. The exact
        // currently positioned graph is merged by the linear backend too.
        const double pitch = found->totalArcLength / found->totalArcElements;
        const double epsilon = std::max(1e-10, pitch * 1e-9);
        std::vector<double> samples{begin, std::nextafter(end, begin)};
        const double firstEvent = std::floor(begin / pitch) * pitch;
        for (double event = firstEvent; event <= end + pitch; event += pitch) {
            if (event >= begin && event < end) {
                samples.push_back(std::max(begin, event - epsilon));
                samples.push_back(std::min(std::nextafter(end, begin), event + epsilon));
            }
        }
        for (double sample = begin; sample < end; sample += pitch * 0.5)
            samples.push_back(sample);
        for (double sample : samples) {
            auto candidate = *found;
            positionGap(candidate, sample, 0.0);
            appendAgeStructure(candidate, definition.upperEntries);
        }
    }
    closeStructureUnderPeriodicConstraints(*m_solver, definition.upperEntries);
    std::sort(definition.upperEntries.begin(), definition.upperEntries.end());
    definition.upperEntries.erase(
        std::unique(definition.upperEntries.begin(), definition.upperEntries.end()),
        definition.upperEntries.end());
    return definition;
}
#endif

void FSolverAnalysisBackend::synchronize(const ModelDefinition &model,
                                         const SolveParameters &parameters,
                                         const PreparedAnalysis &prepared,
                                         std::shared_ptr<const mesh::SolverMesh> mesh,
                                         std::uint64_t topologyIdentity, Dirty rebuilt)
{
    if (!mesh)
        throw std::invalid_argument("FSolver backend requires a mesh");
    configure(model, parameters, prepared);
    const bool topologyChanged = topologyIdentity != m_topologyIdentity;
    if (topologyChanged) {
        if (m_solver->LoadMesh(*mesh) != NOERROR)
            throw std::runtime_error("FSolver could not import the session mesh");
        std::vector<std::pair<std::size_t, std::size_t>> connectivity;
        connectivity.reserve(mesh->edges.size());
        for (const auto &edge : mesh->edges)
            connectivity.emplace_back(edge.first, edge.second);
        if (!m_solver->Cuthill(connectivity))
            throw std::runtime_error("FSolver Cuthill-McKee ordering failed");
        m_topologyIdentity = topologyIdentity;
        m_mesh = std::move(mesh);
        m_lastSystem.reset();
        m_trialSolution.clear();
        m_committedSolution.clear();
        m_haveTrialSolution = false;
        m_haveCommittedSolution = false;
#ifdef XFEMM_USE_CUDSS
        m_currentBucket = {};
        m_cudssBuckets.clear();
#endif
        m_pendingBucketConstructionMs = 0;
        ++m_topologyImports;
        ++m_orderings;
    }
    if (topologyChanged || (rebuilt & Dirty::AirGapCoupling) != Dirty::None)
        positionAirGaps(prepared);
    else
        m_lastAirGapUpdateMs = 0;
#ifdef XFEMM_USE_CUDSS
    if (m_useCudss) {
        const std::string identity = cudssBucketIdentity(prepared);
        auto found = m_cudssBuckets.find(identity);
        if (found == m_cudssBuckets.end()) {
            const auto started = std::chrono::steady_clock::now();
            CudssBucketDefinition definition = buildCudssBucket(prepared);
            found = m_cudssBuckets.emplace(identity, std::move(definition)).first;
            m_pendingBucketConstructionMs +=
                std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - started).count();
        }
        m_currentBucket = found->second;
    }
#endif
}

TrialSolution FSolverAnalysisBackend::solve(const ModelDefinition &model,
                                            const SolveParameters &parameters,
                                            const PreparedAnalysis &prepared)
{
    return solveConfigured(model, parameters, prepared, false);
}

TrialSolution FSolverAnalysisBackend::solveConfigured(
    const ModelDefinition &model, const SolveParameters &parameters,
    const PreparedAnalysis &prepared, bool initialization)
{
    const auto evaluateStarted = std::chrono::steady_clock::now();
    // Static2D/StaticAxisymmetric currently assemble both parts of a fresh
    // matrix on every nonlinear iteration. The persistent cuDSS path retains
    // its bucket graph, analysis, allocations, factor workspace and vectors.
    ++m_operatorAssemblies;
    ++m_rightHandSideAssemblies;
    const auto modelPreparationStarted = std::chrono::steady_clock::now();
    configure(model, parameters, prepared);
    const double modelPreparationMs =
        std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - modelPreparationStarted).count();
    if (parameters.frequency != 0)
        throw std::invalid_argument("FSolverAnalysisBackend currently returns real (zero-frequency) solutions only");

    bool legacyFallback = false;
    const auto createSystem = [&]() {
        if (m_useCudss) {
#ifdef XFEMM_USE_CUDSS
            return std::unique_ptr<femm::LinearSystemBackend<double>>(
                new CudssLinearSystemBackend({m_cudssOptions.deterministic}));
#else
            throw std::runtime_error("cuDSS backend is not available in this build");
#endif
        }
        auto system = femm::create_backend<double>(femm::default_backend_kind());
        if (!system)
            throw std::runtime_error("FSolver could not create the linear system backend");
        return system;
    };

    if (!m_useCudss || !m_lastSystem
#ifdef XFEMM_USE_CUDSS
        || dynamic_cast<CudssLinearSystemBackend *>(m_lastSystem.get()) == nullptr
#endif
        ) {
        m_lastSystem = createSystem();
        if (!m_lastSystem)
            throw std::runtime_error("FSolver could not create the linear system backend");
        m_lastSystem->set_precision(m_solver->Precision);
        if (!m_lastSystem->create(m_solver->NumNodes, m_solver->BandWidth))
            throw std::runtime_error("FSolver could not allocate the linear system");
    }

    femm::LinearSystemBackend<double> &system = *m_lastSystem;
    system.set_precision(m_solver->Precision);
    system.reset_diagnostics();
    const double bucketDefinitionConstructionMs = m_pendingBucketConstructionMs;
    m_pendingBucketConstructionMs = 0;
#ifdef XFEMM_USE_CUDSS
    if (auto *cudss = dynamic_cast<CudssLinearSystemBackend *>(&system))
        cudss->activateBucket(m_currentBucket);
#endif
    if (m_useCudss)
        system.wipe();

    const double outputScale = 4.0 * PI * 1.e-5;
    const std::vector<double> *seed = nullptr;
    if (m_haveTrialSolution) seed = &m_trialSolution;
    else if (m_haveCommittedSolution) seed = &m_committedSolution;
    if (!seed && parameters.initialState && parameters.initialState->realNodalState) {
        const auto &portable = *parameters.initialState->realNodalState;
        if (portable.size() == static_cast<std::size_t>(m_solver->NumNodes)) {
            m_trialSolution.resize(portable.size());
            std::transform(portable.begin(), portable.end(), m_trialSolution.begin(),
                           [outputScale](double value) { return value / outputScale; });
            seed = &m_trialSolution;
        }
    }
    const bool warmStart = m_useCudss && seed &&
                           seed->size() == static_cast<std::size_t>(m_solver->NumNodes);
    if (warmStart)
        std::copy(seed->begin(), seed->end(), system.solution().begin());
    m_solver->setSessionWarmStartUsed(warmStart);

    bool solved = false;
    try {
        solved = m_solver->ProblemType == PLANAR
               ? m_solver->Static2D(system) : m_solver->StaticAxisymmetric(system);
    } catch (const std::exception &error) {
        if (!m_useCudss || !m_cudssOptions.explicitLegacyFallback)
            throw std::runtime_error(std::string("persistent cuDSS evaluation failed: ") +
                                     error.what());
        std::fprintf(stderr,
                     "xfemm: persistent cuDSS evaluation failed; explicitly restarting "
                     "this evaluation with legacy CPU backend: %s\n", error.what());
        legacyFallback = true;
        m_lastSystem = femm::create_backend<double>(BackendKind::Legacy);
        if (!m_lastSystem ||
            !m_lastSystem->create(m_solver->NumNodes, m_solver->BandWidth))
            throw std::runtime_error("legacy fallback allocation failed");
        m_lastSystem->set_precision(m_solver->Precision);
        if (warmStart)
            std::copy(seed->begin(), seed->end(), m_lastSystem->solution().begin());
        m_solver->setSessionWarmStartUsed(warmStart);
        solved = m_solver->ProblemType == PLANAR
               ? m_solver->Static2D(*m_lastSystem)
               : m_solver->StaticAxisymmetric(*m_lastSystem);
    }
    if (!solved)
        throw std::runtime_error(m_useCudss
            ? "persistent cuDSS linear/nonlinear solve did not converge"
            : "FSolver failed to solve the analysis");
    ++m_solves;

    femm::LinearSystemBackend<double> &solvedSystem = *m_lastSystem;
    m_trialSolution.assign(solvedSystem.solution().begin(),
                           solvedSystem.solution().end());
    m_haveTrialSolution = true;

    const auto resultPackagingStarted = std::chrono::steady_clock::now();
    TrialSolution result;
    result.real.emplace();
    result.real->nodal.magneticVectorPotential.reserve(m_solver->NumNodes);
    for (int i = 0; i < m_solver->NumNodes; ++i)
        result.real->nodal.magneticVectorPotential.push_back(solvedSystem.rhs()[i]);
    result.real->nodal.x.reserve(m_solver->NumNodes);
    result.real->nodal.y.reserve(m_solver->NumNodes);
    for (const auto &node : m_solver->meshnode) {
        // FSolver stores imported coordinates internally in centimetres.
        result.real->nodal.x.push_back(node.x / 100.0);
        result.real->nodal.y.push_back(node.y / 100.0);
    }
    for (std::size_t i = 0;
         i < static_cast<std::size_t>(m_solver->NumCircPropsOrig); ++i) {
        const auto &constraint = parameters.circuitConstraints.at(CircuitId{i});
        CComplex current = constraint.kind == CircuitConstraintKind::PrescribedCurrent
                         ? constraint.value : m_solver->circproplist[i].Amps;
        std::optional<CComplex> voltage;
        if (constraint.kind == CircuitConstraintKind::PrescribedVoltage)
            voltage = constraint.value;
        result.circuits.push_back({CircuitId{i}, current, CComplex(), voltage});
        std::optional<double> realVoltage;
        if (voltage)
            realVoltage = voltage->re;
        result.real->circuits.push_back(
            {CircuitId{i}, current.re, 0.0, realVoltage,
             std::nullopt, std::nullopt});
    }
    const double resultPackagingMs =
        std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - resultPackagingStarted).count();

    const auto linear = solvedSystem.diagnostics();
    const auto staticTimings = m_solver->staticSolveTimings();
    result.diagnostics.sessionInitializationMs = m_initializationMs;
    result.diagnostics.modelPreparationMs = modelPreparationMs;
    result.diagnostics.bucketLookupMs = linear.bucketLookupMs;
    result.diagnostics.bucketSwitchMs = linear.bucketSwitchMs;
    result.diagnostics.bucketDefinitionConstructionMs =
        bucketDefinitionConstructionMs;
    result.diagnostics.bucketConstructionMs =
        bucketDefinitionConstructionMs + linear.bucketConstructionMs;
    result.diagnostics.symbolicAnalysisMs = linear.symbolicAnalysisMs;
    result.diagnostics.airGapUpdateMs = m_lastAirGapUpdateMs;
    result.diagnostics.nonlinearMaterialEvaluationMs =
        staticTimings.nonlinearMaterialEvaluationMs;
    result.diagnostics.numericMatrixAssemblyMs =
        staticTimings.numericMatrixAssemblyMs;
    result.diagnostics.sparsePackingMs = linear.sparsePackingMs;
    result.diagnostics.hostToDeviceMs = linear.hostToDeviceMs;
    result.diagnostics.numericFactorizationMs = linear.numericFactorizationMs;
    result.diagnostics.linearSolveMs = linear.solveMs;
    result.diagnostics.deviceToHostMs = linear.deviceToHostMs;
    result.diagnostics.residualEvaluationMs = linear.residualEvaluationMs;
    result.diagnostics.nonlinearBookkeepingMs =
        staticTimings.nonlinearBookkeepingMs;
    result.diagnostics.resultPackagingMs = resultPackagingMs;
    result.diagnostics.matrixNonzeros = linear.matrixNonzeros;
    result.diagnostics.factorNonzeros = linear.factorNonzeros;
    result.diagnostics.permanentDeviceBytes = linear.permanentDeviceBytes;
    result.diagnostics.peakDeviceBytes = linear.peakDeviceBytes;
    result.diagnostics.hostToDeviceBytes = linear.hostToDeviceBytes;
    result.diagnostics.nonlinearIterations = m_solver->newtonIterations();
    result.diagnostics.relativeResidual = linear.lastRelativeResidual;
    result.diagnostics.linearSolver = linear.solver.empty()
        ? (legacyFallback ? "legacy-pcg-explicit-fallback" : "legacy-pcg")
        : linear.solver;
    result.diagnostics.bucketIdentity = linear.bucketIdentity;
    result.diagnostics.converged = legacyFallback ? true : linear.allConverged;
    result.diagnostics.symbolicReused = linear.symbolicReused;
    result.diagnostics.bucketReused = linear.bucketReused;
    result.diagnostics.nonlinearWarmStartUsed = warmStart;
    result.diagnostics.factorizationRetained = linear.factorizationRetained;
    result.diagnostics.exactTopologyFallback = linear.exactTopologyFallback;
    result.diagnostics.legacyFallback = legacyFallback;
    result.diagnostics.deterministic = linear.deterministic;
    result.diagnostics.totalEvaluateMs = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - evaluateStarted).count();
    if (initialization)
        result.backendState = "cudss-session-initialization";
    return result;
}

void FSolverAnalysisBackend::initialize(const ModelDefinition &model,
                                        const SolveParameters &parameters,
                                        const PreparedAnalysis &prepared)
{
    if (!m_useCudss || m_initialized)
        return;
    const auto started = std::chrono::steady_clock::now();
    TrialSolution initialized = solveConfigured(model, parameters, prepared, true);
    m_committedSolution = m_trialSolution;
    m_haveCommittedSolution = m_haveTrialSolution;
    m_initializationMs += std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - started).count();
    m_initializationDiagnostics = initialized.diagnostics;
    m_initializationDiagnostics.sessionInitializationMs = m_initializationMs;
    m_initialized = true;
}

void FSolverAnalysisBackend::commitTrial(const TrialSolution &)
{
    if (!m_useCudss || !m_haveTrialSolution)
        return;
    m_committedSolution = m_trialSolution;
    m_haveCommittedSolution = true;
}

void FSolverAnalysisBackend::rollbackToCommitted()
{
    if (!m_useCudss)
        return;
    if (!m_haveCommittedSolution)
        throw std::logic_error("there is no committed magnetic state to restore");
    m_trialSolution = m_committedSolution;
    m_haveTrialSolution = true;
    if (m_lastSystem &&
        m_lastSystem->dimension() == static_cast<int>(m_trialSolution.size())) {
        std::copy(m_trialSolution.begin(), m_trialSolution.end(),
                  m_lastSystem->solution().begin());
        const double outputScale = 4.0 * PI * 1.e-5;
        for (std::size_t i = 0; i < m_trialSolution.size(); ++i)
            m_lastSystem->rhs()[i] = m_trialSolution[i] * outputScale;
    }
}

void FSolverAnalysisBackend::writeSolution(const std::string &ansPath)
{
    if (!m_lastSystem)
        throw std::logic_error("there is no solved field to export");
    if (ansPath.size() < 4 || ansPath.substr(ansPath.size() - 4) != ".ans")
        throw std::invalid_argument("solution export path must end in .ans");
    m_solver->PathName = ansPath.substr(0, ansPath.size() - 4);
    // The legacy writer handles both planar and axisymmetric static fields.
    const bool written = m_solver->WriteStatic2D(*m_lastSystem);
    if (!written)
        throw std::runtime_error("FSolver could not export the session solution");
}

const FSolver &FSolverAnalysisBackend::solvedSolver() const
{
    if (!m_lastSystem) throw std::logic_error("there is no solved field; call solve first");
    return *m_solver;
}

const femm::LinearSystemBackend<double> &FSolverAnalysisBackend::solvedSystem() const
{
    if (!m_lastSystem) throw std::logic_error("there is no solved field; call solve first");
    return *m_lastSystem;
}

} // namespace femm
