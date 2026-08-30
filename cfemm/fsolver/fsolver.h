/*
   This code is a modified version of an algorithm
   forming part of the software program Finite
   Element Method Magnetics (FEMM), authored by
   David Meeker. The original software code is
   subject to the Aladdin Free Public Licence
   version 8, November 18, 1999. For more information
   on FEMM see www.femm.info. This modified version
   is not endorsed in any way by the original
   authors of FEMM.

   This software has been modified to use the C++
   standard template libraries and remove all Microsoft (TM)
   MFC dependent code to allow easier reuse across
   multiple operating system platforms.

   Date Modified: 2017
   By: Richard Crozier
       Johannes Zarl-Zierl
   Contact:
        richard.crozier@yahoo.co.uk
        johannes.zarl-zierl@jku.at

   Contributions by Johannes Zarl-Zierl were funded by
   Linz Center of Mechatronics GmbH (LCM)
*/

// fsolver.h : interface of the FSolver class
//
/////////////////////////////////////////////////////////////////////////////

#ifndef FSOLVER_H
#define FSOLVER_H

#include <array>
#include <string>
#include <vector>
#include "feasolver.h"
#include "linsolve/LinearSystemBackend.h"
#include "CBlockLabel.h"
#include "CCircuit.h"
#include "CElement.h"
#include "CMaterialProp.h"
#include "CNode.h"
#include "CPointProp.h"
#include "mesh/SolverMesh.h"

namespace femm {
class LuaInstance;
}

class FSolver : public FEASolver<
        femm::CMPointProp
        , femm::CMBoundaryProp
        , femm::CMSolverMaterialProp
        , femm::CMCircuit
        , femm::CMBlockLabel
        , femmsolver::CMElement
        >
{

// Attributes
public:

    struct StaticSolveTimings
    {
        double nonlinearMaterialEvaluationMs = 0;
        double numericMatrixAssemblyMs = 0;
        double airGapMatrixAssemblyMs = 0;
        double elementMatrixAssemblyMs = 0;
        double rhsConstructionMs = 0;
        double boundaryConditionApplicationMs = 0;
        double nonlinearBookkeepingMs = 0;
    };

    struct SweepWarmStartState
    {
        std::vector<femm::CNode> nodes;
        std::vector<double> solution;
    };

    FSolver();
    ~FSolver();

    // General problem attributes
    double Frequency;  ///< \brief Frequency for harmonic problems [Hz]
    double  Relax;

    // mesh information
    std::vector <femm::CNode> meshnode;
    int NumCircPropsOrig;


// Operations
public:

    LoadMeshErr LoadMesh(bool deleteFiles=true) override;
    /** Import a session-owned, value-based mesh without serializing it to disk. */
    LoadMeshErr LoadMesh(const femm::mesh::SolverMesh &mesh);
    /**
     * @brief loadPreviousSolution
     * @return \c true on success, \c false otherwise.
     * \internal
     * ### FEMM reference source
     *  - \femm42{fkn/femmedoccore.cpp,CFemmeDocCore::LoadPrev()}
     * \endinternal
     */
    bool loadPreviousSolution(bool loadAprev);
    //bool LoadMeshFromPrevSolution(bool loadAprev);
    bool LoadMeshNodesFromSolution(bool loadA, FILE* fp);
    bool LoadMeshElementsFromSolution(FILE* fp);
    bool LoadPBCFromSolution(FILE* fp);
    bool LoadAGEsFromSolution(FILE* fp);
    bool LoadProblemFile();
    int Static2D(femm::LinearSystemBackend<double> &L);
    /**
     * @brief WriteStatic2D
     * @param L
     * @return \c true on success, \c false otherwise.
     * \internal
     * ### FEMM reference source
     *  - \femm42{fkn/prob1big.cpp,CFemmeDocCore::WriteStatic2D()}
     * \endinternal
     */
    int WriteStatic2D(femm::LinearSystemBackend<double> &L);
    int Harmonic2D(femm::LinearSystemBackend<CComplex> &L,bool verbose=false);
    int WriteHarmonic2D(femm::LinearSystemBackend<CComplex> &L);
    int StaticAxisymmetric(femm::LinearSystemBackend<double> &L);
    int HarmonicAxisymmetric(femm::LinearSystemBackend<CComplex> &L,bool verbose=false);
    void GetFillFactor(int lbl);
    double ElmArea(int i);

    virtual bool runSolver(bool verbose=false) override;

    void setSweepWarmStart(const SweepWarmStartState &state)
    {
        sweepInitialState = state;
    }
    void setPreparedMesh(const femm::mesh::SolverMesh &mesh)
    {
        preparedMesh = &mesh;
    }
    const SweepWarmStartState &acceptedSweepState() const
    {
        return sweepAcceptedState;
    }
    bool usedSweepWarmStart() const { return sweepWarmStartUsed; }
    bool remappedSweepWarmStart() const { return sweepWarmStartRemapped; }
    int newtonIterations() const { return lastNewtonIterations; }
    double finalNewtonRelativeUpdate() const { return lastNewtonRelativeUpdate; }
    const StaticSolveTimings &staticSolveTimings() const { return lastStaticSolveTimings; }
    /** Apply the homogeneous part of Static2D's Dirichlet and periodic
     * constraints to column-major sensitivity right-hand sides. */
    void transformStatic2DSensitivityRhs(std::vector<double> &rightHandSides,
                                         std::size_t rightHandSideCount) const;
    /** Session-only seed marker; direct solvers ignore it as a linear initial guess. */
    void setSessionWarmStartUsed(bool used) { sweepWarmStartUsed = used; }

private:

    struct Static2DElementAssemblyData
    {
        std::array<int, 3> nodes{};
        std::array<double, 3> p{};
        std::array<double, 3> q{};
        std::array<double, 9> mx{};
        std::array<double, 9> my{};
        std::array<double, 9> mxy{};
        std::array<double, 9> fixedMatrix{};
        std::array<double, 3> fixedRhs{};
        double area = 0;
        int block = -1;
        int label = -1;
        int circuit = -1;
    };

    struct Static2DDynamicAssemblyData
    {
        std::array<double, 9> nonlinearMatrix{};
        std::array<double, 6> upperMatrix{};
        std::array<double, 3> rhs{};
    };

    void prepareStatic2DAssemblyData();
    void invalidateStatic2DAssemblyData()
    {
        static2DElementAssemblyData.clear();
        static2DDynamicAssemblyData.clear();
        static2DPlan = {};
    }

    virtual void CleanUp() override;

    /**
     * @brief getPrevAxiB
     * @param k
     * @param B1p
     * @param B2p
     * \internal
     * ### FEMM reference source
     *  - \femm42{fkn/prob4big.cpp,CFemmeDocCore::GetPrevAxiB()}
     * \endinternal
     */
    void getPrevAxiB(int k, double &B1p, double &B2p) const;
    /**
     * @brief getPrev2DB
     * @param k
     * @param B1p
     * @param B2p
     * \internal
     * ### FEMM reference source
     *  - \femm42{fkn/prob2big.cpp,CFemmeDocCore::GetPrev2DB()}
     * \endinternal
     */
    void getPrev2DB(int k, double &B1p, double &B2p) const;

    // override parent class virtual method
    void SortNodes (std::vector<int> newnum) override;

    bool handleToken(const std::string &token, std::istream &input, std::ostream &err) override;

    femm::LuaInstance *theLua;

    /// Vector containing previous solution for incremental permeability analysis
    std::vector <double> Aprev;

    SweepWarmStartState sweepInitialState;
    SweepWarmStartState sweepAcceptedState;
    std::vector<double> sweepMappedInitialSolution;
    const femm::mesh::SolverMesh *preparedMesh = nullptr;
    bool sweepWarmStartUsed = false;
    bool sweepWarmStartRemapped = false;
    int lastNewtonIterations = 0;
    double lastNewtonRelativeUpdate = 0;
    StaticSolveTimings lastStaticSolveTimings;
    std::vector<Static2DElementAssemblyData> static2DElementAssemblyData;
    std::vector<Static2DDynamicAssemblyData> static2DDynamicAssemblyData;
    femm::PlanarAssemblyPlan static2DPlan;
};

/////////////////////////////////////////////////////////////////////////////

double GetNewMu(double mu,int BHpoints, CComplex *BHdata,double muc,double B);

#endif
