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

//		CFknDlg *TheView;

private:
    /**
     * Pack the linked-list assembly representation into contiguous upper
     * triangular rows.  Assembly benefits from stable entry addresses, while
     * PCG traverses the matrix thousands of times and must not pointer-chase
     * individually allocated entries on every pass.
     */
    void rebuildCompactRows();
    bool solveParallelPCG(int flag);
    void applyParallelPreconditionerChunk(const double *X, double *Y,
                                          int begin, int end);

    std::vector<std::size_t> m_rowOffsets;
    std::vector<int> m_columns;
    std::vector<double> m_values;
    std::vector<std::size_t> m_symmetricRowOffsets;
    std::vector<int> m_symmetricColumns;
    std::vector<double> m_symmetricValues;
    bool m_compactRowsDirty = true;
    int m_parallelThreads = 1;
    bool m_useJacobi = false;
    bool m_useParallelPcg = false;
};

#endif
