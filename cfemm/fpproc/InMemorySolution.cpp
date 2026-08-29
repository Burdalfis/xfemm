#include "fpproc.h"

#include "fsolver.h"
#include "linsolve/LinearSystemBackend.h"

#include <stdexcept>
#include <cmath>
#include <cstdlib>

namespace {

std::vector<int> copyPhysicalLabels(const FSolver &solver,
                                    std::vector<femm::CMBlockLabel> &target)
{
    std::vector<int> mapping(solver.labellist.size(), -1);
    target.clear();
    target.reserve(solver.labellist.size());
    for (std::size_t i = 0; i < solver.labellist.size(); ++i) {
        if (solver.labellist[i].BlockType < 0) continue;
        mapping[i] = static_cast<int>(target.size());
        target.push_back(solver.labellist[i]);
        auto &label = target.back();
        if (label.InCircuit >= 0 &&
            label.InCircuit < static_cast<int>(solver.circproplist.size())) {
            const int original = solver.circproplist[
                static_cast<std::size_t>(label.InCircuit)].OrigCirc;
            if (original >= 0)
                label.InCircuit = original;
        }
    }
    return mapping;
}

void copyOriginalCircuits(const FSolver &solver,
                          std::vector<femm::CMCircuit> &target,
                          const femm::FemmProblem *problem)
{
    const std::size_t count = static_cast<std::size_t>(solver.NumCircPropsOrig);
    if (count > solver.circproplist.size())
        throw std::runtime_error("solver original circuit count is invalid");
    std::vector<int> circuitTypes;
    circuitTypes.reserve(count);
    if (problem) {
        if (problem->circproplist.size() != count)
            throw std::runtime_error("model/solver circuit count mismatch");
        for (const auto &property : problem->circproplist) {
            const auto *circuit = dynamic_cast<const femm::CMCircuit *>(property.get());
            if (!circuit)
                throw std::runtime_error("magnetic model contains a non-magnetic circuit");
            circuitTypes.push_back(circuit->CircType);
        }
    } else {
        if (target.size() != count)
            throw std::runtime_error("persistent postprocessor circuit count changed");
        for (const auto &circuit : target)
            circuitTypes.push_back(circuit.CircType);
    }
    target.assign(solver.circproplist.begin(),
                  solver.circproplist.begin() + static_cast<std::ptrdiff_t>(count));
    for (std::size_t i = 0; i < count; ++i)
        target[i].CircType = circuitTypes[i];
}

void copyPhysicalCircuitState(const FSolver &solver,
                              const std::vector<int> &mapping,
                              std::vector<femm::CMBlockLabel> &target)
{
    for (std::size_t source = 0; source < solver.labellist.size(); ++source) {
        if (source >= mapping.size() || mapping[source] < 0)
            continue;
        auto &label = target[static_cast<std::size_t>(mapping[source])];
        const int circuit = solver.labellist[source].InCircuit;
        if (circuit < 0) {
            label.Case = 1;
            label.J = 0.;
        } else {
            if (circuit >= static_cast<int>(solver.circproplist.size()))
                throw std::runtime_error("solver block label circuit is invalid");
            const auto &solvedCircuit = solver.circproplist[
                static_cast<std::size_t>(circuit)];
            if (solvedCircuit.Case == 0) {
                label.Case = 0;
                label.dVolts = solvedCircuit.dVolts;
            } else {
                label.Case = 1;
                label.J = solvedCircuit.J;
            }
        }
    }
}

void releaseAirGapFields(femmsolver::CAirGapElement &age)
{
    std::free(age.brc); std::free(age.brs);
    std::free(age.btc); std::free(age.bts);
    std::free(age.br); std::free(age.bt);
    std::free(age.brcPrev); std::free(age.brsPrev);
    std::free(age.btcPrev); std::free(age.btsPrev);
    std::free(age.brPrev); std::free(age.btPrev);
    std::free(age.nh);
    age.brc = age.brs = age.btc = age.bts = nullptr;
    age.br = age.bt = nullptr;
    age.brcPrev = age.brsPrev = age.btcPrev = age.btsPrev = nullptr;
    age.brPrev = age.btPrev = nullptr;
    age.nh = nullptr;
}

void clearAirGapFieldPointers(femmsolver::CAirGapElement &age)
{
    age.brc = age.brs = age.btc = age.bts = nullptr;
    age.br = age.bt = nullptr;
    age.brcPrev = age.brsPrev = age.btcPrev = age.btsPrev = nullptr;
    age.brPrev = age.btPrev = nullptr;
    age.nh = nullptr;
}

} // namespace

bool FPProc::OpenDocument(const femm::FemmProblem &problem, const FSolver &solver,
                          const femm::LinearSystemBackend<double> &solution)
{
    if (solution.dimension() != solver.NumNodes)
        throw std::invalid_argument("solution and solver node counts differ");

    NewDocument();
    Frequency = solver.Frequency;
    Depth = problem.Depth;
    Precision = problem.Precision;
    LengthUnits = problem.LengthUnits;
    problemType = problem.problemType;
    Coords = problem.Coords;
    ProblemNote = problem.comment;
    extRo = problem.extRo;
    extRi = problem.extRi;
    extZo = problem.extZo;
    PrevSoln = problem.previousSolutionFile;
    PrevType = problem.PrevType;
    bIncremental = problem.PrevType;

    for (const auto &item : problem.nodelist) nodelist.push_back(*item);
    for (const auto &item : problem.linelist) linelist.push_back(*item);
    for (const auto &item : problem.arclist) arclist.push_back(*item);
    const auto labelMapping = copyPhysicalLabels(solver, blocklist);
    nodeproplist = solver.nodeproplist;
    lineproplist = solver.lineproplist;
    blockproplist.clear();
    for (const auto &source : solver.blockproplist) {
        blockproplist.push_back(static_cast<const femm::CMMaterialProp &>(source));
        auto &item = blockproplist.back();
        item.MuMax = problem.PrevType ? source.MuMax : 0.;
        item.Frequency = source.Frequency;
        item.mu_fdx = source.mu_fdx;
        item.mu_fdy = source.mu_fdy;
        item.clearSlopes();
        item.GetSlopes(2. * std::acos(-1.) * Frequency);
        if (item.BHpoints > 0 && item.slope.size() < static_cast<std::size_t>(item.BHpoints))
            throw std::runtime_error("could not prepare nonlinear material slopes");
    }
    copyOriginalCircuits(solver, circproplist, &problem);

    const double centimetresPerSourceUnit[] = {2.54, 0.1, 1., 100., 0.00254, 1.e-04};
    const double coordinateScale = centimetresPerSourceUnit[LengthUnits];
    meshnode.resize(solver.meshnode.size());
    for (std::size_t i = 0; i < solver.meshnode.size(); ++i) {
        meshnode[i].x = solver.meshnode[i].x / coordinateScale;
        meshnode[i].y = solver.meshnode[i].y / coordinateScale;
        meshnode[i].A = solution.rhs()[i];
    }

    meshelem.resize(solver.meshele.size());
    for (std::size_t i = 0; i < solver.meshele.size(); ++i) {
        static_cast<femmsolver::CMElement &>(meshelem[i]) = solver.meshele[i];
        if (meshelem[i].lbl < 0 ||
            meshelem[i].lbl >= static_cast<int>(labelMapping.size()) ||
            labelMapping[static_cast<std::size_t>(meshelem[i].lbl)] < 0)
            throw std::runtime_error("solver element refers to a non-physical label");
        meshelem[i].lbl = labelMapping[static_cast<std::size_t>(meshelem[i].lbl)];
        meshelem[i].blk = blocklist[meshelem[i].lbl].BlockType;
    }

    // The solver stores AGE lengths in centimetres; FPProc uses metres here.
    agelist = solver.agelist;
    for (auto &age : agelist) {
        age.ri /= 100.;
        age.ro /= 100.;
        age.agc /= 100.;
    }
    NumAirGapElems = static_cast<int>(agelist.size());
    pmeshnode = &meshnode;
    pmeshelem = &meshelem;

    copyPhysicalCircuitState(solver, labelMapping, blocklist);

    return finalizeSolution();
}

bool FPProc::UpdateSolution(const FSolver &solver,
                            const femm::LinearSystemBackend<double> &solution)
{
    if (solution.dimension() != solver.NumNodes ||
        meshnode.size() != solver.meshnode.size() ||
        meshelem.size() != solver.meshele.size())
        throw std::invalid_argument("persistent postprocessor topology changed");

    Frequency = solver.Frequency;
    Precision = solver.Precision;
    const auto labelMapping = copyPhysicalLabels(solver, blocklist);
    copyOriginalCircuits(solver, circproplist, nullptr);
    for (std::size_t i = 0; i < meshnode.size(); ++i)
        meshnode[i].A = solution.rhs()[i];
    for (std::size_t i = 0; i < meshelem.size(); ++i) {
        static_cast<femmsolver::CMElement &>(meshelem[i]) = solver.meshele[i];
        if (meshelem[i].lbl < 0 ||
            meshelem[i].lbl >= static_cast<int>(labelMapping.size()) ||
            labelMapping[static_cast<std::size_t>(meshelem[i].lbl)] < 0)
            throw std::runtime_error("solver element refers to a non-physical label");
        meshelem[i].lbl = labelMapping[static_cast<std::size_t>(meshelem[i].lbl)];
        meshelem[i].blk = blocklist[meshelem[i].lbl].BlockType;
    }

    for (auto &age : agelist) releaseAirGapFields(age);
    agelist = solver.agelist;
    for (auto &age : agelist) {
        clearAirGapFieldPointers(age);
        age.ri /= 100.;
        age.ro /= 100.;
        age.agc /= 100.;
    }
    NumAirGapElems = static_cast<int>(agelist.size());

    copyPhysicalCircuitState(solver, labelMapping, blocklist);
    return finalizeSolution(false);
}

bool FPProc::UpdateAirGapSolution(
    const FSolver &solver,
    const femm::LinearSystemBackend<double> &solution)
{
    if (solution.dimension() != solver.NumNodes ||
        meshnode.size() != solver.meshnode.size())
        throw std::invalid_argument("persistent postprocessor topology changed");

    Frequency = solver.Frequency;
    Precision = solver.Precision;
    for (std::size_t i = 0; i < meshnode.size(); ++i)
        meshnode[i].A = solution.rhs()[i];

    for (auto &age : agelist) releaseAirGapFields(age);
    agelist = solver.agelist;
    for (auto &age : agelist) {
        clearAirGapFieldPointers(age);
        age.ri /= 100.;
        age.ro /= 100.;
        age.agc /= 100.;
    }
    NumAirGapElems = static_cast<int>(agelist.size());
    return finalizeSolution(false, false);
}
