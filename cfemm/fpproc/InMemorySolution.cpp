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
    }
    return mapping;
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
    circproplist = solver.circproplist;

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

    for (std::size_t i = 0; i < blocklist.size(); ++i) {
        const int circuit = blocklist[i].InCircuit;
        if (circuit < 0) {
            blocklist[i].Case = 1;
            blocklist[i].J = 0.;
        } else if (circproplist[circuit].Case == 0) {
            blocklist[i].Case = 0;
            blocklist[i].dVolts = circproplist[circuit].dVolts;
        } else {
            blocklist[i].Case = 1;
            blocklist[i].J = circproplist[circuit].J;
        }
    }

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
    circproplist = solver.circproplist;
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

    for (std::size_t i = 0; i < blocklist.size(); ++i) {
        const int circuit = blocklist[i].InCircuit;
        if (circuit < 0) {
            blocklist[i].Case = 1;
            blocklist[i].J = 0.;
        } else if (circproplist[circuit].Case == 0) {
            blocklist[i].Case = 0;
            blocklist[i].dVolts = circproplist[circuit].dVolts;
        } else {
            blocklist[i].Case = 1;
            blocklist[i].J = circproplist[circuit].J;
        }
    }
    return finalizeSolution(false);
}
