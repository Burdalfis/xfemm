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

   Date Modified: 2011 - 11 - 10
   By: Richard Crozier
   Contact: richard.crozier@yahoo.co.uk
*/

#ifndef SPARS_H
#define SPARS_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

class CEntry
{
public:

    double x;				// value stored in the entry
    int c;					// column that the entry lives in
    CEntry *next;			// pointer to next entry in row;
    CEntry();

private:
};


class CBigLinProb
{
public:

    // data members

    double *V;				// solution
    double *P;				// search direction;
    double *R;				// residual;
    double *U;				// A * P;
    double *Z;
    double *b;				// RHS of linear equation
    CEntry **M;				// pointer to list of matrix entries;
    int n;					// dimensions of the matrix;
    int bdw;				// Optional matrix bandwidth parameter;
    double Precision;		// error tolerance for solution
    double Lambda;			// relaxation factor;

    int *Q; ///< Used by esolver and hsolver.

    // member functions

    // constructor
    CBigLinProb();
    // destructor
    ~CBigLinProb();
    virtual int Create(int d, int bw);	// initialize the problem
    void Put(double v, int p, int q);
    // use to create/set entries in the matrix
    double Get(int p, int q);
    bool PCGSolve(int flag);	// flag==true if guess for V present;
    void MultPC(const double *X, double *Y);
    void AddTo(double v, int p, int q);
    void MultA(const double *X, double *Y);
    void SetValue(int i, double x);
    void Periodicity(int i, int j);
    void AntiPeriodicity(int i, int j);
    void Wipe();
    double Dot(const double *X, const double *Y);
    void ComputeBandwidth();
    int lastIterations() const { return m_lastIterations; }
    double lastRelativeResidual() const { return m_lastRelativeResidual; }

    /**
     * Copy the assembled symmetric matrix's stored upper triangle into
     * conventional zero-based CSR.  This is intentionally a read-only bridge
     * from FEMM's linked-list assembly representation to native sparse
     * libraries.  Explicitly assembled zero entries are retained because
     * they are part of a persistent solver's immutable sparsity graph.
     */
    void copyUpperCsr(std::vector<std::int32_t> &rowOffsets,
                      std::vector<std::int32_t> &columnIndices,
                      std::vector<double> &values);

//		CFknDlg *TheView;

private:
    /**
     * Pack the linked-list assembly representation into contiguous upper
     * triangular rows.  Assembly benefits from stable entry addresses, while
     * PCG traverses the matrix thousands of times and must not pointer-chase
     * individually allocated entries on every pass.
     */
    void rebuildCompactRows();
    void writeLinearSystemExport(const std::string &prefix,
                                 unsigned long long exportIndex,
                                 bool warmStart,
                                 const std::vector<double> &initialSolution,
                                 bool byTopology);
    bool solveParallelPCG(int flag);
    void applyParallelPreconditionerChunk(const double *X, double *Y,
                                          int begin, int end);
    std::size_t mixedWideRowOffset(int row) const
    {
        return m_wideRowOffsets16.empty()
            ? static_cast<std::size_t>(m_wideRowOffsets32[static_cast<std::size_t>(row)])
            : static_cast<std::size_t>(m_wideRowOffsets16[static_cast<std::size_t>(row)]);
    }

    std::vector<std::size_t> m_rowOffsets;
    std::vector<int> m_columns;
    std::vector<std::uint16_t> m_columnOffsets16;
    std::vector<std::uint16_t> m_wideRowOffsets16;
    std::vector<std::uint32_t> m_wideRowOffsets32;
    std::vector<int> m_wideColumns;
    std::vector<int> m_wideRowIndices;
    std::vector<double> m_values;
    std::vector<std::size_t> m_symmetricRowOffsets;
    std::vector<int> m_symmetricColumns;
    std::vector<double> m_symmetricValues;
    bool m_compactRowsDirty = true;
    int m_parallelThreads = 1;
    bool m_useJacobi = false;
    bool m_useParallelPcg = false;
    bool m_useMixedColumnOffsets = false;
    bool m_useRowColumnOffsets = false;
    bool m_collectStats = false;
    bool m_matrixStatsPrinted = false;
    int m_lastIterations = 0;
    double m_lastRelativeResidual = -1.;
    unsigned long long m_solveCount = 0;
    unsigned long long m_totalIterations = 0;
    unsigned long long m_spmvCalls = 0;
    unsigned long long m_preconditionerCalls = 0;
    unsigned long long m_dotCalls = 0;
    unsigned long long m_packCalls = 0;
    double m_spmvSeconds = 0.;
    double m_preconditionerSeconds = 0.;
    double m_dotSeconds = 0.;
    double m_packSeconds = 0.;
};

#endif
