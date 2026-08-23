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

#include "femmcomplex.h"
#include "spars.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <utility>

#ifdef XFEMM_USE_OPENMP
#include <omp.h>
#endif

using std::swap;

#define KLUDGE


CEntry::CEntry()
{
    next=NULL;
    x=0;
    c=0;
}

CBigLinProb::CBigLinProb()
{
    n=0;
    // Best guess for relaxation parameter
    Lambda = 1.5;

#ifdef XFEMM_USE_OPENMP
    if (const char *threads = std::getenv("XFEMM_NUM_THREADS"))
    {
        char *end = nullptr;
        const long requested = std::strtol(threads, &end, 10);
        if (end != threads && *end == '\0' && requested > 0 && requested <= 1024)
            m_parallelThreads = static_cast<int>(requested);
    }
    // Contiguous block-SSOR is the threaded default. It preserves most of
    // SSOR's local coupling without introducing cross-thread dependencies.
    m_useParallelPcg = m_parallelThreads > 1;
    if (const char *preconditioner = std::getenv("XFEMM_PCG_PRECONDITIONER"))
    {
        if (std::strcmp(preconditioner, "jacobi") == 0)
        {
            m_useJacobi = true;
            m_useParallelPcg = m_parallelThreads > 1;
        }
        else if (std::strcmp(preconditioner, "block-ssor") == 0)
        {
            m_useJacobi = false;
            m_useParallelPcg = m_parallelThreads > 1;
        }
        else if (std::strcmp(preconditioner, "ssor") == 0)
        {
            m_useJacobi = false;
            m_useParallelPcg = false;
        }
    }
#endif
}

CBigLinProb::~CBigLinProb()
{
    if (n==0) return;

    int i;
    CEntry *uo,*ui;

    free(b);
    free(P);
    free(R);
    free(V);
    free(U);
    free(Z);

    for(i=0; i<n; i++)
    {
        ui=M[i];
        do
        {
            uo=ui;
            ui=uo->next;
            delete uo;
        }
        while(ui!=NULL);
    }

    free(M);
    free(Q);
    n = 0;
}

int CBigLinProb::Create(int d, int bw)
{
    int i;

    bdw=bw;
    b=(double *)calloc(d,sizeof(double));
    V=(double *)calloc(d,sizeof(double));
    P=(double *)calloc(d,sizeof(double));
    R=(double *)calloc(d,sizeof(double));
    U=(double *)calloc(d,sizeof(double));
    Z=(double *)calloc(d,sizeof(double));

    M=(CEntry **)calloc(d,sizeof(CEntry *));
    n=d;

    for(i=0; i<d; i++)
    {
        M[i] = new CEntry;
        M[i]->c = i;
    }
    Q = (int *)  calloc(d,sizeof(int));

    return 1;
}

void CBigLinProb::Put(double v, int p, int q)
{
    CEntry *e,*l = NULL;

    m_compactRowsDirty = true;

    if (q<p)
        swap(p,q);

    e = M[p];

    while ((e->c < q) && (e->next != NULL))
    {
        l = e;
        e = e->next;
    }

    if (e->c == q)
    {
        e->x = v;
        return;
    }

    CEntry *m = new CEntry;

    if ((e->next == NULL) && (q > e->c))
    {
        e->next = m;
        m->c = q;
        m->x = v;
    }
    else
    {
        l->next = m;
        m->next = e;
        m->c = q;
        m->x = v;
    }
    return;
}

double CBigLinProb::Get(int p, int q)
{
    if (q < p)
    {
        swap(p,q);
    }

    CEntry *e = M[p];
    while ((e->c < q) && (e->next != NULL))
    {
        e = e->next;
    }

    if (e->c == q) return e->x;

    return 0;
}

void CBigLinProb::AddTo(double v, int p, int q)
{
	Put(Get(p,q)+v,p,q);
}

void CBigLinProb::rebuildCompactRows()
{
    if (!m_compactRowsDirty)
        return;

    m_rowOffsets.resize(static_cast<std::size_t>(n) + 1);
    m_columns.clear();
    m_values.clear();

    for (int i = 0; i < n; ++i)
    {
        m_rowOffsets[static_cast<std::size_t>(i)] = m_values.size();
        for (CEntry *entry = M[i]; entry != NULL; entry = entry->next)
        {
            m_columns.push_back(entry->c);
            m_values.push_back(entry->x);
        }
    }
    m_rowOffsets[static_cast<std::size_t>(n)] = m_values.size();

    // Only the parallel path needs both halves. Scalar SSOR/SpMV uses the
    // more compact upper triangle and avoids the extra memory traffic.
    if (!m_useParallelPcg)
    {
        m_symmetricRowOffsets.clear();
        m_symmetricColumns.clear();
        m_symmetricValues.clear();
        m_compactRowsDirty = false;
        return;
    }

    // Materialize both halves for race-free row-oriented SpMV. Processing
    // source rows in ascending order also appends every destination row in
    // ascending column order, matching the legacy accumulation order.
    m_symmetricRowOffsets.assign(static_cast<std::size_t>(n) + 1, 0);
    for (int i = 0; i < n; ++i)
    {
        const std::size_t begin = m_rowOffsets[static_cast<std::size_t>(i)];
        const std::size_t end = m_rowOffsets[static_cast<std::size_t>(i) + 1];
        for (std::size_t j = begin; j < end; ++j)
        {
            ++m_symmetricRowOffsets[static_cast<std::size_t>(i) + 1];
            if (m_columns[j] != i)
                ++m_symmetricRowOffsets[static_cast<std::size_t>(m_columns[j]) + 1];
        }
    }
    for (int i = 0; i < n; ++i)
        m_symmetricRowOffsets[static_cast<std::size_t>(i) + 1] +=
            m_symmetricRowOffsets[static_cast<std::size_t>(i)];

    m_symmetricColumns.resize(m_symmetricRowOffsets.back());
    m_symmetricValues.resize(m_symmetricRowOffsets.back());
    std::vector<std::size_t> next = m_symmetricRowOffsets;
    for (int i = 0; i < n; ++i)
    {
        const std::size_t begin = m_rowOffsets[static_cast<std::size_t>(i)];
        const std::size_t end = m_rowOffsets[static_cast<std::size_t>(i) + 1];
        for (std::size_t j = begin; j < end; ++j)
        {
            const int column = m_columns[j];
            std::size_t destination = next[static_cast<std::size_t>(i)]++;
            m_symmetricColumns[destination] = column;
            m_symmetricValues[destination] = m_values[j];
            if (column != i)
            {
                destination = next[static_cast<std::size_t>(column)]++;
                m_symmetricColumns[destination] = i;
                m_symmetricValues[destination] = m_values[j];
            }
        }
    }
    m_compactRowsDirty = false;
}

void CBigLinProb::MultA(const double *X, double *Y)
{
    rebuildCompactRows();

    for(int i=0; i<n; i++) Y[i]=0.;

    for(int i=0; i<n; i++)
    {
        const std::size_t begin = m_rowOffsets[static_cast<std::size_t>(i)];
        const std::size_t end = m_rowOffsets[static_cast<std::size_t>(i) + 1];
        Y[i] += m_values[begin] * X[i];
        for(std::size_t j = begin + 1; j < end; ++j)
        {
            const int column = m_columns[j];
            const double value = m_values[j];
            Y[i] += value * X[column];
            Y[column] += value * X[i];
        }
    }
}

double CBigLinProb::Dot(const double *X, const double *Y)
{
    double z=0.;

    for(int i=0; i<n; i++) z+=X[i]*Y[i];

    return z;
}

void CBigLinProb::MultPC(const double *X, double *Y)
{
    rebuildCompactRows();

    if (m_useJacobi)
    {
        for (int i = 0; i < n; ++i)
            Y[i] = X[i] / m_values[m_rowOffsets[static_cast<std::size_t>(i)]];
        return;
    }

    const double c = Lambda*(2.-Lambda);
    for(int i=0; i<n; i++) Y[i]=X[i]*c;

    // invert Lower Triangle;
    for(int i=0; i<n; i++)
    {
        const std::size_t begin = m_rowOffsets[static_cast<std::size_t>(i)];
        const std::size_t end = m_rowOffsets[static_cast<std::size_t>(i) + 1];
        Y[i] /= m_values[begin];
        for(std::size_t j = begin + 1; j < end; ++j)
            Y[m_columns[j]] -= m_values[j] * Y[i] * Lambda;
    }

    for(int i=0; i<n; i++)
        Y[i] *= m_values[m_rowOffsets[static_cast<std::size_t>(i)]];

    // invert Upper Triangle
    for(int i=n-1; i>=0; i--)
    {
        const std::size_t begin = m_rowOffsets[static_cast<std::size_t>(i)];
        const std::size_t end = m_rowOffsets[static_cast<std::size_t>(i) + 1];
        for(std::size_t j = begin + 1; j < end; ++j)
            Y[i] -= m_values[j] * Y[m_columns[j]] * Lambda;
        Y[i] /= m_values[begin];
    }
}

bool CBigLinProb::PCGSolve(int flag)
{
    int i;
    double res,res_o,res_new;
    double er,del,rho,pAp;

    // quick check for most obvious sign of singularity;
    for(i=0; i<n; i++) if(M[i]->x==0)
        {
            fprintf(stderr,"singular flag tripped at %i of %i\n", i,n);
            return 0;
        }

    // initialize progress bar;
//	TheView->SetDlgItemText(IDC_FRAME1,"Conjugate Gradient Solver");
//	TheView->m_prog1.SetPos(0);
    printf("Conjugate Gradient Solver\n");

#ifdef XFEMM_USE_OPENMP
    if (m_useParallelPcg)
        return solveParallelPCG(flag);
#endif

    // residual with V=0
    MultPC(b,Z);
    res_o=Dot(Z,b);
    if(res_o==0) return true;

    // if flag is false, initialize V with zeros;
    if (flag==0) for(i=0; i<n; i++) V[i]=0;

    // form residual;
    MultA(V,R);
    for(i=0; i<n; i++) R[i]=b[i]-R[i];

    // form initial search direction;
    MultPC(R,Z);
    for(i=0; i<n; i++) P[i]=Z[i];
    res=Dot(Z,R);

    // do iteration;
    do
    {
        // step i)
        MultA(P,U);
        pAp=Dot(P,U);
        del=res/pAp;

        for(i=0; i<n; i++)
        {
            // step ii)
            V[i]+=(del*P[i]);

            // step iii)
            R[i]-=(del*U[i]);
        }

        // step iv)
        MultPC(R,Z);
        res_new=Dot(Z,R);
        rho=res_new/res;
        res=res_new;

        // step v)
        for(i=0; i<n; i++) P[i]=Z[i]+(rho*P[i]);

        // have we converged yet?
        er=sqrt(res/res_o);
//        prg2=(int) (20.*log10(er)/(log10(Precision)));
//        if(prg2>prg1)
//        {
//            prg1=prg2;
//            prg2=(prg1*5);
//            if(prg2>100) prg2=100;
//			TheView->m_prog1.SetPos(prg2);
//			TheView->InvalidateRect(NULL, FALSE);
//			TheView->UpdateWindow();
//        }

    }
    while(er>Precision);

    return true;
}

#ifdef XFEMM_USE_OPENMP
void CBigLinProb::applyParallelPreconditionerChunk(const double *X, double *Y,
                                                   int begin, int end)
{
    if (m_useJacobi)
    {
        for (int i = begin; i < end; ++i)
            Y[i] = X[i] / m_values[m_rowOffsets[static_cast<std::size_t>(i)]];
        return;
    }

    const double scale = Lambda * (2. - Lambda);
    for (int i = begin; i < end; ++i)
        Y[i] = X[i] * scale;

    for (int i = begin; i < end; ++i)
    {
        const std::size_t rowBegin = m_rowOffsets[static_cast<std::size_t>(i)];
        const std::size_t rowEnd = m_rowOffsets[static_cast<std::size_t>(i) + 1];
        Y[i] /= m_values[rowBegin];
        for (std::size_t j = rowBegin + 1; j < rowEnd; ++j)
        {
            const int column = m_columns[j];
            if (column >= end)
                break;
            Y[column] -= m_values[j] * Y[i] * Lambda;
        }
    }

    for (int i = begin; i < end; ++i)
        Y[i] *= m_values[m_rowOffsets[static_cast<std::size_t>(i)]];

    for (int i = end - 1; i >= begin; --i)
    {
        const std::size_t rowBegin = m_rowOffsets[static_cast<std::size_t>(i)];
        const std::size_t rowEnd = m_rowOffsets[static_cast<std::size_t>(i) + 1];
        for (std::size_t j = rowBegin + 1; j < rowEnd; ++j)
        {
            const int column = m_columns[j];
            if (column >= end)
                break;
            Y[i] -= m_values[j] * Y[column] * Lambda;
        }
        Y[i] /= m_values[rowBegin];
    }
}

bool CBigLinProb::solveParallelPCG(int flag)
{
    rebuildCompactRows();

    double res_o = 0.;
    double res = 0.;
    double res_new = 0.;
    double pAp = 0.;
    double del = 0.;
    double rho = 0.;
    double er = 0.;
    bool active = true;

#pragma omp parallel num_threads(m_parallelThreads) \
    shared(res_o, res, res_new, pAp, del, rho, er, active)
    {
        const int thread = omp_get_thread_num();
        const int threadCount = omp_get_num_threads();
        const int blockBegin = static_cast<int>(
            (static_cast<long long>(n) * thread) / threadCount);
        const int blockEnd = static_cast<int>(
            (static_cast<long long>(n) * (thread + 1)) / threadCount);

        applyParallelPreconditionerChunk(b, Z, blockBegin, blockEnd);
#pragma omp barrier
#pragma omp for schedule(static) reduction(+:res_o)
        for (int i = 0; i < n; ++i)
            res_o += Z[i] * b[i];

#pragma omp single
        active = res_o != 0.;

        if (active)
        {
            if (flag == 0)
            {
#pragma omp for schedule(static)
                for (int i = 0; i < n; ++i)
                    V[i] = 0.;
            }

#pragma omp for schedule(static)
            for (int i = 0; i < n; ++i)
            {
                const std::size_t begin =
                    m_symmetricRowOffsets[static_cast<std::size_t>(i)];
                const std::size_t end =
                    m_symmetricRowOffsets[static_cast<std::size_t>(i) + 1];
                double sum = 0.;
                for (std::size_t j = begin; j < end; ++j)
                    sum += m_symmetricValues[j] * V[m_symmetricColumns[j]];
                R[i] = b[i] - sum;
            }

            applyParallelPreconditionerChunk(R, Z, blockBegin, blockEnd);
#pragma omp barrier
#pragma omp single
            res = 0.;
#pragma omp for schedule(static) reduction(+:res)
            for (int i = 0; i < n; ++i)
            {
                P[i] = Z[i];
                res += Z[i] * R[i];
            }

            do
            {
#pragma omp single
                pAp = 0.;
#pragma omp for schedule(static) reduction(+:pAp)
                for (int i = 0; i < n; ++i)
                {
                    const std::size_t begin =
                        m_symmetricRowOffsets[static_cast<std::size_t>(i)];
                    const std::size_t end =
                        m_symmetricRowOffsets[static_cast<std::size_t>(i) + 1];
                    double sum = 0.;
                    for (std::size_t j = begin; j < end; ++j)
                        sum += m_symmetricValues[j] * P[m_symmetricColumns[j]];
                    U[i] = sum;
                    pAp += P[i] * sum;
                }

#pragma omp single
                del = res / pAp;

#pragma omp for schedule(static)
                for (int i = 0; i < n; ++i)
                {
                    V[i] += del * P[i];
                    R[i] -= del * U[i];
                }

                applyParallelPreconditionerChunk(R, Z, blockBegin, blockEnd);
#pragma omp barrier
#pragma omp single
                res_new = 0.;
#pragma omp for schedule(static) reduction(+:res_new)
                for (int i = 0; i < n; ++i)
                    res_new += Z[i] * R[i];

#pragma omp single
                {
                    rho = res_new / res;
                    res = res_new;
                    er = sqrt(res / res_o);
                }

#pragma omp for schedule(static)
                for (int i = 0; i < n; ++i)
                    P[i] = Z[i] + rho * P[i];
            }
            while (er > Precision);
        }
    }

    return true;
}
#endif

void CBigLinProb::SetValue(int i, double x)
{
    int k,fst,lst;
    double z;

    if(bdw==0)
    {
        fst=0;
        lst=n;
    }
    else
    {
        fst=i-bdw;
        if (fst<0) fst=0;
        lst=i+bdw;
        if (lst>n) lst=n;
    }

    for(k=fst; k<lst; k++)
    {
        z=Get(k,i);
        if(z!=0)
        {
            b[k]=b[k]-(z*x);
            if(i!=k) Put(0.,k,i);
        }
    }
    b[i]=Get(i,i)*x;
}

void CBigLinProb::Wipe()
{
    int i;
    CEntry *e;

    for(i=0; i<n; i++)
    {
        b[i]=0.;
        e=M[i];
        do
        {
            e->x=0;
            e=e->next;
        }
        while(e!=NULL);
    }
    m_compactRowsDirty = true;
}

void CBigLinProb::AntiPeriodicity(int i, int j)
{
    int k,fst,lst;
    double v1,v2,c;

#ifdef KLUDGE
    int tmpbdw=bdw;
    bdw=0;
#endif

    if (j<i)
        swap(j,i);

    if(bdw==0)
    {
        fst=0;
        lst=n;
    }
    else
    {
        fst=i-bdw;
        if (fst<0) fst=0;
        lst=j+bdw;
        if (lst>n) lst=n;
    }

    for(k=fst; k<lst; k++)
    {
        if((k!=i) && (k!=j))
        {
            v1=Get(k,i);
            v2=Get(k,j);
            if ((v1!=0) || (v2!=0))
            {
                c=(v1-v2)/2.;
                Put(c,k,i);
                Put(-c,k,j);
            }
        }
        if((k==i+bdw) && (k<j-bdw) && (bdw!=0)) k=j-bdw;
    }

    c=0.5*(Get(i,i)+Get(j,j));
    Put(c,i,i);
    Put(c,j,j);

    c=0.5*(b[i]-b[j]);
    b[i]=c;
    b[j]=-c;

#ifdef KLUDGE
    bdw=tmpbdw;
#endif
}

void CBigLinProb::Periodicity(int i, int j)
{
    int k,fst,lst;
    double v1,v2,c;

#ifdef KLUDGE
    int tmpbdw=bdw;
    bdw=0;
#endif

    if (j<i)
        swap(j,i);

    if(bdw==0)
    {
        fst=0;
        lst=n;
    }
    else
    {
        fst=i-bdw;
        if (fst<0) fst=0;
        lst=j+bdw;
        if (lst>n) lst=n;
    }

    for(k=fst; k<lst; k++)
    {
        if((k!=i) && (k!=j))
        {
            v1=Get(k,i);
            v2=Get(k,j);
            if ((v1!=0) || (v2!=0))
            {
                c=(v1+v2)/2.;
                Put(c,k,i);
                Put(c,k,j);
            }
        }
        if((k==i+bdw) && (k<j-bdw) && (bdw!=0)) k=j-bdw;
    }

    c=(Get(i,i)+Get(j,j))/2.;
    Put(c,i,i);
    Put(c,j,j);

    c=0.5*(b[i]+b[j]);
    b[i]=c;
    b[j]=c;

#ifdef KLUDGE
    bdw=tmpbdw;
#endif
}


// a diagnostic routine to check whether that the bandwidth of the
// constructed matrix is actually consistent with a priori bandwidth.
void CBigLinProb::ComputeBandwidth()
{
    CEntry *e;
    int k,bw,maxbw;

    for(maxbw=0,k=0; k<n; k++)
    {
        e=M[k];
        while(e->next != NULL) e=e->next;
        bw=e->c - k;
        if (bw>maxbw) maxbw=bw;
    }

//	MsgBox("Assumed Bandwidth = %i\nActual Bandwidth = %i",bdw,maxbw);

    printf("Assumed Bandwidth = %i\nActual Bandwidth = %i", bdw, maxbw);
}
