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

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#ifdef XFEMM_USE_OPENMP
#include <omp.h>
#endif

using std::swap;

#define KLUDGE

namespace {

std::atomic<unsigned long long> nextLinearSystemExport{0};

template <typename T>
bool writeBinary(FILE *stream, const T &value)
{
    return std::fwrite(&value, sizeof(value), 1, stream) == 1;
}

template <typename T>
bool writeBinaryArray(FILE *stream, const T *values, std::size_t count)
{
    return count == 0 || std::fwrite(values, sizeof(T), count, stream) == count;
}

void hashBytes(std::uint64_t &hash, const void *data, std::size_t bytes)
{
    const unsigned char *input = static_cast<const unsigned char *>(data);
    for (std::size_t i = 0; i < bytes; ++i)
    {
        hash ^= input[i];
        hash *= 1099511628211ull;
    }
}

} // namespace


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

    if (const char *stats = std::getenv("XFEMM_PCG_STATS"))
        m_collectStats = stats[0] != '\0' && std::strcmp(stats, "0") != 0;

    if (const char *relaxation = std::getenv("XFEMM_PCG_LAMBDA"))
    {
        char *end = nullptr;
        const double requested = std::strtod(relaxation, &end);
        if (end != relaxation && *end == '\0' && std::isfinite(requested) &&
            requested > 0. && requested < 2.)
        {
            Lambda = requested;
        }
        else
        {
            std::fprintf(stderr,
                         "Ignoring invalid XFEMM_PCG_LAMBDA='%s'; expected "
                         "a finite value between 0 and 2\n",
                         relaxation);
        }
    }

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
#endif

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

    if (const char *columnIndex = std::getenv("XFEMM_PCG_COLUMN_INDEX"))
    {
        if (std::strcmp(columnIndex, "mixed16") == 0)
            m_useMixedColumnOffsets = true;
        else if (std::strcmp(columnIndex, "row16") == 0)
            m_useRowColumnOffsets = true;
        else if (columnIndex[0] != '\0' && std::strcmp(columnIndex, "int32") != 0)
            std::fprintf(stderr,
                         "Ignoring invalid XFEMM_PCG_COLUMN_INDEX='%s'; "
                         "expected int32, mixed16, or row16\n",
                         columnIndex);
    }

    if ((m_useMixedColumnOffsets || m_useRowColumnOffsets) && m_useParallelPcg)
    {
        std::fprintf(stderr,
                     "16-bit XFEMM_PCG_COLUMN_INDEX modes apply only to "
                     "scalar PCG; using int32 for the parallel solver\n");
        m_useMixedColumnOffsets = false;
        m_useRowColumnOffsets = false;
    }
}

CBigLinProb::~CBigLinProb()
{
    if (m_collectStats && m_solveCount != 0)
    {
        std::fprintf(stderr,
                     "XFEMM_PCG_TOTAL solves=%llu iterations=%llu lambda=%.9g "
                     "pack_s=%.6f spmv_s=%.6f preconditioner_s=%.6f "
                     "dot_s=%.6f\n",
                     m_solveCount, m_totalIterations, Lambda, m_packSeconds,
                     m_spmvSeconds, m_preconditionerSeconds, m_dotSeconds);
    }

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

    const auto packStarted = m_collectStats
        ? std::chrono::steady_clock::now()
        : std::chrono::steady_clock::time_point{};

    m_rowOffsets.resize(static_cast<std::size_t>(n) + 1);
    m_columns.clear();
    m_columnOffsets16.clear();
    m_wideRowOffsets16.clear();
    m_wideRowOffsets32.clear();
    m_wideColumns.clear();
    m_wideRowIndices.clear();
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

    if (m_useMixedColumnOffsets)
    {
        const std::uint16_t escape = std::numeric_limits<std::uint16_t>::max();
        m_columnOffsets16.resize(m_columns.size());
        std::vector<std::uint32_t> wideRowOffsets(
            static_cast<std::size_t>(n) + 1, 0);
        for (int i = 0; i < n; ++i)
        {
            const std::size_t begin = m_rowOffsets[static_cast<std::size_t>(i)];
            const std::size_t end = m_rowOffsets[static_cast<std::size_t>(i) + 1];
            wideRowOffsets[static_cast<std::size_t>(i)] =
                static_cast<std::uint32_t>(m_wideColumns.size());
            for (std::size_t j = begin; j < end; ++j)
            {
                const int offset = m_columns[j] - i;
                if (offset >= 0 && offset < static_cast<int>(escape))
                {
                    m_columnOffsets16[j] = static_cast<std::uint16_t>(offset);
                }
                else
                {
                    m_columnOffsets16[j] = escape;
                    m_wideColumns.push_back(m_columns[j]);
                }
            }
        }
        wideRowOffsets[static_cast<std::size_t>(n)] =
            static_cast<std::uint32_t>(m_wideColumns.size());

        if (m_wideColumns.size() <=
            static_cast<std::size_t>(std::numeric_limits<std::uint16_t>::max()))
        {
            m_wideRowOffsets16.resize(wideRowOffsets.size());
            std::transform(wideRowOffsets.begin(), wideRowOffsets.end(),
                           m_wideRowOffsets16.begin(),
                           [](std::uint32_t offset) {
                               return static_cast<std::uint16_t>(offset);
                           });
        }
        else
        {
            m_wideRowOffsets32.swap(wideRowOffsets);
        }
    }
    else if (m_useRowColumnOffsets)
    {
        m_columnOffsets16.reserve(m_columns.size());
        for (int i = 0; i < n; ++i)
        {
            const std::size_t begin = m_rowOffsets[static_cast<std::size_t>(i)];
            const std::size_t end = m_rowOffsets[static_cast<std::size_t>(i) + 1];
            bool wideRow = false;
            for (std::size_t j = begin + 1; j < end; ++j)
            {
                if (m_columns[j] - i >
                    static_cast<int>(std::numeric_limits<std::uint16_t>::max()))
                {
                    wideRow = true;
                    break;
                }
            }
            if (wideRow)
                m_wideRowIndices.push_back(i);
            for (std::size_t j = begin; j < end; ++j)
            {
                if (wideRow)
                    m_wideColumns.push_back(m_columns[j]);
                else
                    m_columnOffsets16.push_back(
                        static_cast<std::uint16_t>(m_columns[j] - i));
            }
        }
    }

    if (m_collectStats && !m_matrixStatsPrinted)
    {
        std::vector<int> bandwidths;
        bandwidths.reserve(m_columns.size() > static_cast<std::size_t>(n)
                               ? m_columns.size() - static_cast<std::size_t>(n)
                               : 0);
        unsigned long long bandwidthSum = 0;
        std::size_t wideBandwidthCount = 0;
        std::size_t wideRowCount = 0;
        std::size_t wideRowNnz = 0;
        int maximumBandwidth = 0;
        for (int i = 0; i < n; ++i)
        {
            const std::size_t begin = m_rowOffsets[static_cast<std::size_t>(i)];
            const std::size_t end = m_rowOffsets[static_cast<std::size_t>(i) + 1];
            bool wideRow = false;
            for (std::size_t j = begin + 1; j < end; ++j)
            {
                const int bandwidth = m_columns[j] - i;
                bandwidths.push_back(bandwidth);
                bandwidthSum += static_cast<unsigned long long>(bandwidth);
                if (bandwidth > 65535)
                {
                    ++wideBandwidthCount;
                    wideRow = true;
                }
                maximumBandwidth = std::max(maximumBandwidth, bandwidth);
            }
            if (wideRow)
            {
                ++wideRowCount;
                wideRowNnz += end - begin;
            }
        }
        std::sort(bandwidths.begin(), bandwidths.end());
        const auto percentile = [&bandwidths](double fraction) {
            if (bandwidths.empty())
                return 0;
            const std::size_t index = static_cast<std::size_t>(
                fraction * static_cast<double>(bandwidths.size() - 1));
            return bandwidths[index];
        };
        std::size_t columnBytes = m_columns.size() * sizeof(m_columns[0]);
        if (m_useRowColumnOffsets)
        {
            columnBytes =
                m_columnOffsets16.size() * sizeof(m_columnOffsets16[0]) +
                m_wideColumns.size() * sizeof(m_wideColumns[0]) +
                m_wideRowIndices.size() * sizeof(m_wideRowIndices[0]);
        }
        else if (m_useMixedColumnOffsets)
        {
            columnBytes =
                m_columnOffsets16.size() * sizeof(m_columnOffsets16[0]) +
                m_wideRowOffsets16.size() * sizeof(m_wideRowOffsets16[0]) +
                m_wideRowOffsets32.size() * sizeof(m_wideRowOffsets32[0]) +
                m_wideColumns.size() * sizeof(m_wideColumns[0]);
        }
        const std::size_t packedBytes =
            m_rowOffsets.size() * sizeof(m_rowOffsets[0]) +
            columnBytes + m_values.size() * sizeof(m_values[0]);
        const std::size_t rowHybridBytesWithoutTags =
            m_rowOffsets.size() * sizeof(m_rowOffsets[0]) +
            m_values.size() * sizeof(m_values[0]) +
            (m_columns.size() - wideRowNnz) * sizeof(std::uint16_t) +
            wideRowNnz * sizeof(m_columns[0]);
        const std::size_t rowHybridTaggedBytes =
            rowHybridBytesWithoutTags + static_cast<std::size_t>(n);
        const std::size_t rowHybridBytes =
            rowHybridBytesWithoutTags + wideRowCount * sizeof(int);
        const double meanBandwidth = bandwidths.empty()
            ? 0.
            : static_cast<double>(bandwidthSum) /
                  static_cast<double>(bandwidths.size());
        std::fprintf(
            stderr,
            "XFEMM_PCG_MATRIX n=%d upper_nnz=%zu offdiag_nnz=%zu lambda=%.9g "
            "packed_bytes=%zu bandwidth_max=%d bandwidth_mean=%.3f "
            "bandwidth_p50=%d bandwidth_p90=%d bandwidth_p95=%d "
            "bandwidth_p99=%d bandwidth_over_uint16=%zu "
            "bandwidth_over_uint16_pct=%.6f uint16_offsets=%s "
            "wide_rows=%zu wide_rows_pct=%.6f wide_row_nnz=%zu "
            "wide_row_nnz_pct=%.6f row_hybrid_bytes_no_row_map=%zu "
            "row_hybrid_tagged_bytes=%zu row_hybrid_bytes=%zu "
            "column_index=%s wide_columns=%zu "
            "wide_row_offset_bits=%d\n",
            n, m_values.size(), bandwidths.size(), Lambda, packedBytes,
            maximumBandwidth, meanBandwidth, percentile(0.50),
            percentile(0.90), percentile(0.95), percentile(0.99),
            wideBandwidthCount,
            bandwidths.empty()
                ? 0.
                : 100. * static_cast<double>(wideBandwidthCount) /
                      static_cast<double>(bandwidths.size()),
            maximumBandwidth <= 65535 ? "yes" : "no",
            wideRowCount,
            n == 0 ? 0. : 100. * static_cast<double>(wideRowCount) /
                                  static_cast<double>(n),
            wideRowNnz,
            m_columns.empty()
                ? 0.
                : 100. * static_cast<double>(wideRowNnz) /
                      static_cast<double>(m_columns.size()),
            rowHybridBytesWithoutTags, rowHybridTaggedBytes, rowHybridBytes,
            m_useRowColumnOffsets
                ? "row16"
                : (m_useMixedColumnOffsets ? "mixed16" : "int32"),
            m_wideColumns.size(),
            m_useMixedColumnOffsets
                ? (m_wideRowOffsets16.empty() ? 32 : 16)
                : 0);
        m_matrixStatsPrinted = true;
    }

    // Only the parallel path needs both halves. Scalar SSOR/SpMV uses the
    // more compact upper triangle and avoids the extra memory traffic.
    if (!m_useParallelPcg)
    {
        m_symmetricRowOffsets.clear();
        m_symmetricColumns.clear();
        m_symmetricValues.clear();
        if (m_useMixedColumnOffsets || m_useRowColumnOffsets)
            std::vector<int>().swap(m_columns);
        m_compactRowsDirty = false;
        if (m_collectStats)
        {
            ++m_packCalls;
            m_packSeconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - packStarted).count();
        }
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
    if (m_collectStats)
    {
        ++m_packCalls;
        m_packSeconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - packStarted).count();
    }
}

void CBigLinProb::copyUpperCsr(std::vector<std::int32_t> &rowOffsets,
                               std::vector<std::int32_t> &columnIndices,
                               std::vector<double> &values)
{
    rebuildCompactRows();
    if (m_rowOffsets.back() >
        static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()))
        throw std::overflow_error("real sparse matrix has more than INT32_MAX stored entries");

    rowOffsets.resize(m_rowOffsets.size());
    std::transform(m_rowOffsets.begin(), m_rowOffsets.end(), rowOffsets.begin(),
                   [](std::size_t offset) {
                       return static_cast<std::int32_t>(offset);
                   });
    columnIndices.assign(m_columns.begin(), m_columns.end());
    values = m_values;
}

void CBigLinProb::writeLinearSystemExport(
    const std::string &prefix, unsigned long long exportIndex, bool warmStart,
    const std::vector<double> &initialSolution, bool byTopology)
{
    // The development format is deliberately simple and fixed-width.  It is
    // documented alongside the standalone CUDA benchmark reader.
    rebuildCompactRows();

    const std::vector<std::size_t> *rowOffsets = &m_symmetricRowOffsets;
    const std::vector<int> *columns = &m_symmetricColumns;
    const std::vector<double> *values = &m_symmetricValues;
    std::vector<std::size_t> fullRowOffsets;
    std::vector<int> fullColumns;
    std::vector<double> fullValues;

    if (rowOffsets->empty())
    {
        fullRowOffsets.assign(static_cast<std::size_t>(n) + 1, 0);
        for (int row = 0; row < n; ++row)
        {
            for (CEntry *entry = M[row]; entry != nullptr; entry = entry->next)
            {
                ++fullRowOffsets[static_cast<std::size_t>(row) + 1];
                if (entry->c != row)
                    ++fullRowOffsets[static_cast<std::size_t>(entry->c) + 1];
            }
        }
        for (int row = 0; row < n; ++row)
            fullRowOffsets[static_cast<std::size_t>(row) + 1] +=
                fullRowOffsets[static_cast<std::size_t>(row)];

        fullColumns.resize(fullRowOffsets.back());
        fullValues.resize(fullRowOffsets.back());
        std::vector<std::size_t> next = fullRowOffsets;
        for (int row = 0; row < n; ++row)
        {
            for (CEntry *entry = M[row]; entry != nullptr; entry = entry->next)
            {
                const int column = entry->c;
                std::size_t destination = next[static_cast<std::size_t>(row)]++;
                fullColumns[destination] = column;
                fullValues[destination] = entry->x;
                if (column != row)
                {
                    destination = next[static_cast<std::size_t>(column)]++;
                    fullColumns[destination] = row;
                    fullValues[destination] = entry->x;
                }
            }
        }
        rowOffsets = &fullRowOffsets;
        columns = &fullColumns;
        values = &fullValues;
    }

    std::vector<std::uint64_t> offsets64(rowOffsets->size());
    std::transform(rowOffsets->begin(), rowOffsets->end(), offsets64.begin(),
                   [](std::size_t offset) {
                       return static_cast<std::uint64_t>(offset);
                   });
    std::uint64_t topologyHash = 1469598103934665603ull;
    hashBytes(topologyHash, offsets64.data(),
              offsets64.size() * sizeof(offsets64[0]));
    hashBytes(topologyHash, columns->data(),
              columns->size() * sizeof((*columns)[0]));

    char suffix[64];
    if (byTopology)
        std::snprintf(suffix, sizeof(suffix), ".%016llx.xfemm-system",
                      static_cast<unsigned long long>(topologyHash));
    else
        std::snprintf(suffix, sizeof(suffix), ".%06llu.xfemm-system",
                      exportIndex);
    const std::string path = prefix + suffix;
    // The indexed temporary name also keeps concurrent same-topology exports
    // from sharing an intermediate file; rename still publishes atomically.
    const std::string temporaryPath =
        path + "." + std::to_string(exportIndex) + ".tmp";
    FILE *stream = std::fopen(temporaryPath.c_str(), "wb");
    if (stream == nullptr)
    {
        std::fprintf(stderr, "XFEMM_SYSTEM_EXPORT error=open path=%s\n",
                     temporaryPath.c_str());
        return;
    }

    const char magic[8] = {'X', 'F', 'E', 'M', 'M', 'L', 'S', '\0'};
    const std::uint32_t version = 1;
    const std::uint32_t endianMarker = 0x01020304u;
    const std::uint32_t scalarBytes = sizeof(double);
    const std::uint32_t columnIndexBytes = sizeof(std::int32_t);
    const std::uint64_t dimension = static_cast<std::uint64_t>(n);
    const std::uint64_t nonzeros = static_cast<std::uint64_t>(values->size());
    const std::uint64_t solveIndex = exportIndex;
    const std::uint32_t flags = (warmStart ? 1u : 0u) | 2u; // full symmetric CSR
    const std::uint32_t parallelThreads =
        static_cast<std::uint32_t>(m_parallelThreads);
    const std::int64_t iterations = static_cast<std::int64_t>(m_lastIterations);

    bool ok = writeBinaryArray(stream, magic, sizeof(magic)) &&
              writeBinary(stream, version) && writeBinary(stream, endianMarker) &&
              writeBinary(stream, scalarBytes) &&
              writeBinary(stream, columnIndexBytes) &&
              writeBinary(stream, dimension) && writeBinary(stream, nonzeros) &&
              writeBinary(stream, solveIndex) && writeBinary(stream, flags) &&
              writeBinary(stream, parallelThreads) &&
              writeBinary(stream, Precision) && writeBinary(stream, Lambda) &&
              writeBinary(stream, iterations) &&
              writeBinary(stream, m_lastRelativeResidual);

    ok = ok && writeBinaryArray(stream, offsets64.data(), offsets64.size());
    ok = ok && writeBinaryArray(stream, columns->data(), columns->size());
    ok = ok && writeBinaryArray(stream, values->data(), values->size());
    ok = ok && writeBinaryArray(stream, b, static_cast<std::size_t>(n));
    ok = ok && writeBinaryArray(stream, initialSolution.data(),
                                initialSolution.size());
    ok = ok && writeBinaryArray(stream, V, static_cast<std::size_t>(n));
    const bool closeOk = std::fclose(stream) == 0;
    ok = ok && closeOk;

    if (!ok || std::rename(temporaryPath.c_str(), path.c_str()) != 0)
    {
        std::remove(temporaryPath.c_str());
        std::fprintf(stderr, "XFEMM_SYSTEM_EXPORT error=write path=%s\n",
                     path.c_str());
        return;
    }
    std::fprintf(stderr,
                 "XFEMM_SYSTEM_EXPORT path=%s n=%d nnz=%zu solve=%llu "
                 "warm_start=%s iterations=%d relative_residual=%.9g "
                 "topology_hash=%016llx\n",
                 path.c_str(), n, values->size(), exportIndex,
                 warmStart ? "yes" : "no", m_lastIterations,
                 m_lastRelativeResidual,
                 static_cast<unsigned long long>(topologyHash));
}

void CBigLinProb::MultA(const double *X, double *Y)
{
    rebuildCompactRows();
    const auto started = m_collectStats
        ? std::chrono::steady_clock::now()
        : std::chrono::steady_clock::time_point{};

    for(int i=0; i<n; i++) Y[i]=0.;

    if (m_useRowColumnOffsets)
    {
        std::size_t narrowIndex = 0;
        std::size_t wideIndex = 0;
        const auto multiplyNarrowRow = [&](int i)
        {
            const std::size_t begin = m_rowOffsets[static_cast<std::size_t>(i)];
            const std::size_t end = m_rowOffsets[static_cast<std::size_t>(i) + 1];
            const std::size_t rowNnz = end - begin;
            Y[i] += m_values[begin] * X[i];
            const std::uint16_t *offsets =
                m_columnOffsets16.data() + narrowIndex;
            for (std::size_t k = 1; k < rowNnz; ++k)
            {
                const int column = i + static_cast<int>(offsets[k]);
                const double value = m_values[begin + k];
                Y[i] += value * X[column];
                Y[column] += value * X[i];
            }
            narrowIndex += rowNnz;
        };
        const auto multiplyWideRow = [&](int i)
        {
            const std::size_t begin = m_rowOffsets[static_cast<std::size_t>(i)];
            const std::size_t end = m_rowOffsets[static_cast<std::size_t>(i) + 1];
            const std::size_t rowNnz = end - begin;
            Y[i] += m_values[begin] * X[i];
            const int *columns = m_wideColumns.data() + wideIndex;
            for (std::size_t k = 1; k < rowNnz; ++k)
            {
                const int column = columns[k];
                const double value = m_values[begin + k];
                Y[i] += value * X[column];
                Y[column] += value * X[i];
            }
            wideIndex += rowNnz;
        };

        int row = 0;
        for (const int wideRow : m_wideRowIndices)
        {
            for (; row < wideRow; ++row)
                multiplyNarrowRow(row);
            multiplyWideRow(row);
            ++row;
        }
        for (; row < n; ++row)
            multiplyNarrowRow(row);
    }
    else if (m_useMixedColumnOffsets)
    {
        const std::uint16_t escape = std::numeric_limits<std::uint16_t>::max();
        for (int i = 0; i < n; ++i)
        {
            const std::size_t begin = m_rowOffsets[static_cast<std::size_t>(i)];
            const std::size_t end = m_rowOffsets[static_cast<std::size_t>(i) + 1];
            std::size_t wide = mixedWideRowOffset(i);
            Y[i] += m_values[begin] * X[i];
            for (std::size_t j = begin + 1; j < end; ++j)
            {
                const std::uint16_t offset = m_columnOffsets16[j];
                const int column = offset == escape
                    ? m_wideColumns[wide++]
                    : i + static_cast<int>(offset);
                const double value = m_values[j];
                Y[i] += value * X[column];
                Y[column] += value * X[i];
            }
        }
    }
    else
    {
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
    if (m_collectStats)
    {
        ++m_spmvCalls;
        m_spmvSeconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count();
    }
}

double CBigLinProb::Dot(const double *X, const double *Y)
{
    const auto started = m_collectStats
        ? std::chrono::steady_clock::now()
        : std::chrono::steady_clock::time_point{};
    double z=0.;

    for(int i=0; i<n; i++) z+=X[i]*Y[i];

    if (m_collectStats)
    {
        ++m_dotCalls;
        m_dotSeconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count();
    }
    return z;
}

void CBigLinProb::MultPC(const double *X, double *Y)
{
    rebuildCompactRows();
    const auto started = m_collectStats
        ? std::chrono::steady_clock::now()
        : std::chrono::steady_clock::time_point{};

    if (m_useJacobi)
    {
        for (int i = 0; i < n; ++i)
            Y[i] = X[i] / m_values[m_rowOffsets[static_cast<std::size_t>(i)]];
        if (m_collectStats)
        {
            ++m_preconditionerCalls;
            m_preconditionerSeconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - started).count();
        }
        return;
    }

    const double c = Lambda*(2.-Lambda);
    for(int i=0; i<n; i++) Y[i]=X[i]*c;

    if (m_useRowColumnOffsets)
    {
        std::size_t narrowIndex = 0;
        std::size_t wideIndex = 0;
        // invert Lower Triangle;
        const auto lowerNarrowRow = [&](int i)
        {
            const std::size_t begin = m_rowOffsets[static_cast<std::size_t>(i)];
            const std::size_t end = m_rowOffsets[static_cast<std::size_t>(i) + 1];
            const std::size_t rowNnz = end - begin;
            Y[i] /= m_values[begin];
            const std::uint16_t *offsets =
                m_columnOffsets16.data() + narrowIndex;
            for (std::size_t k = 1; k < rowNnz; ++k)
            {
                const int column = i + static_cast<int>(offsets[k]);
                Y[column] -= m_values[begin + k] * Y[i] * Lambda;
            }
            narrowIndex += rowNnz;
        };
        const auto lowerWideRow = [&](int i)
        {
            const std::size_t begin = m_rowOffsets[static_cast<std::size_t>(i)];
            const std::size_t end = m_rowOffsets[static_cast<std::size_t>(i) + 1];
            const std::size_t rowNnz = end - begin;
            Y[i] /= m_values[begin];
            const int *columns = m_wideColumns.data() + wideIndex;
            for (std::size_t k = 1; k < rowNnz; ++k)
                Y[columns[k]] -= m_values[begin + k] * Y[i] * Lambda;
            wideIndex += rowNnz;
        };

        int row = 0;
        for (const int wideRow : m_wideRowIndices)
        {
            for (; row < wideRow; ++row)
                lowerNarrowRow(row);
            lowerWideRow(row);
            ++row;
        }
        for (; row < n; ++row)
            lowerNarrowRow(row);

        for (int i = 0; i < n; ++i)
            Y[i] *= m_values[m_rowOffsets[static_cast<std::size_t>(i)]];

        // invert Upper Triangle;
        narrowIndex = m_columnOffsets16.size();
        wideIndex = m_wideColumns.size();
        const auto upperNarrowRow = [&](int i)
        {
            const std::size_t begin = m_rowOffsets[static_cast<std::size_t>(i)];
            const std::size_t end = m_rowOffsets[static_cast<std::size_t>(i) + 1];
            const std::size_t rowNnz = end - begin;
            narrowIndex -= rowNnz;
            const std::uint16_t *offsets =
                m_columnOffsets16.data() + narrowIndex;
            for (std::size_t k = 1; k < rowNnz; ++k)
            {
                const int column = i + static_cast<int>(offsets[k]);
                Y[i] -= m_values[begin + k] * Y[column] * Lambda;
            }
            Y[i] /= m_values[begin];
        };
        const auto upperWideRow = [&](int i)
        {
            const std::size_t begin = m_rowOffsets[static_cast<std::size_t>(i)];
            const std::size_t end = m_rowOffsets[static_cast<std::size_t>(i) + 1];
            const std::size_t rowNnz = end - begin;
            wideIndex -= rowNnz;
            const int *columns = m_wideColumns.data() + wideIndex;
            for (std::size_t k = 1; k < rowNnz; ++k)
                Y[i] -= m_values[begin + k] * Y[columns[k]] * Lambda;
            Y[i] /= m_values[begin];
        };

        row = n - 1;
        for (auto it = m_wideRowIndices.rbegin();
             it != m_wideRowIndices.rend(); ++it)
        {
            for (; row > *it; --row)
                upperNarrowRow(row);
            upperWideRow(row);
            --row;
        }
        for (; row >= 0; --row)
            upperNarrowRow(row);
    }
    else if (m_useMixedColumnOffsets)
    {
        const std::uint16_t escape = std::numeric_limits<std::uint16_t>::max();
        // invert Lower Triangle;
        for (int i = 0; i < n; ++i)
        {
            const std::size_t begin = m_rowOffsets[static_cast<std::size_t>(i)];
            const std::size_t end = m_rowOffsets[static_cast<std::size_t>(i) + 1];
            std::size_t wide = mixedWideRowOffset(i);
            Y[i] /= m_values[begin];
            for (std::size_t j = begin + 1; j < end; ++j)
            {
                const std::uint16_t offset = m_columnOffsets16[j];
                const int column = offset == escape
                    ? m_wideColumns[wide++]
                    : i + static_cast<int>(offset);
                Y[column] -= m_values[j] * Y[i] * Lambda;
            }
        }

        for (int i = 0; i < n; ++i)
            Y[i] *= m_values[m_rowOffsets[static_cast<std::size_t>(i)]];

        // invert Upper Triangle
        for (int i = n - 1; i >= 0; --i)
        {
            const std::size_t begin = m_rowOffsets[static_cast<std::size_t>(i)];
            const std::size_t end = m_rowOffsets[static_cast<std::size_t>(i) + 1];
            std::size_t wide = mixedWideRowOffset(i);
            for (std::size_t j = begin + 1; j < end; ++j)
            {
                const std::uint16_t offset = m_columnOffsets16[j];
                const int column = offset == escape
                    ? m_wideColumns[wide++]
                    : i + static_cast<int>(offset);
                Y[i] -= m_values[j] * Y[column] * Lambda;
            }
            Y[i] /= m_values[begin];
        }
    }
    else
    {
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
    if (m_collectStats)
    {
        ++m_preconditionerCalls;
        m_preconditionerSeconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count();
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

    const auto solveStarted = m_collectStats
        ? std::chrono::steady_clock::now()
        : std::chrono::steady_clock::time_point{};
    const unsigned long long spmvCallsBefore = m_spmvCalls;
    const unsigned long long preconditionerCallsBefore = m_preconditionerCalls;
    const unsigned long long dotCallsBefore = m_dotCalls;
    const double spmvSecondsBefore = m_spmvSeconds;
    const double preconditionerSecondsBefore = m_preconditionerSeconds;
    const double dotSecondsBefore = m_dotSeconds;
    const unsigned long long packCallsBefore = m_packCalls;
    const double packSecondsBefore = m_packSeconds;
    m_lastIterations = 0;
    m_lastRelativeResidual = -1.;

    const char *topologyExportEnvironment =
        std::getenv("XFEMM_LINEAR_SYSTEM_EXPORT_BY_TOPOLOGY");
    const bool exportByTopology = topologyExportEnvironment != nullptr &&
                                  topologyExportEnvironment[0] != '\0';
    const char *exportPrefixEnvironment = exportByTopology
        ? topologyExportEnvironment
        : std::getenv("XFEMM_LINEAR_SYSTEM_EXPORT");
    const bool exportEnabled = exportPrefixEnvironment != nullptr &&
                               exportPrefixEnvironment[0] != '\0';
    const std::string exportPrefix =
        exportEnabled ? exportPrefixEnvironment : std::string();
    const unsigned long long exportIndex = exportEnabled
        ? nextLinearSystemExport.fetch_add(1, std::memory_order_relaxed)
        : 0;
    std::vector<double> initialSolution;
    if (exportEnabled)
        initialSolution.assign(V, V + n);

#ifdef XFEMM_USE_OPENMP
    if (m_useParallelPcg)
    {
        const bool converged = solveParallelPCG(flag);
        if (m_collectStats)
        {
            const double wallSeconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - solveStarted).count();
            ++m_solveCount;
            m_totalIterations += static_cast<unsigned long long>(m_lastIterations);
            std::fprintf(
                stderr,
                "XFEMM_PCG_SOLVE solve=%llu iterations=%d "
                "relative_residual=%.9g wall_s=%.6f parallel=yes "
                "pack_s=%.6f pack_calls=%llu\n",
                m_solveCount - 1, m_lastIterations, m_lastRelativeResidual,
                wallSeconds, m_packSeconds - packSecondsBefore,
                m_packCalls - packCallsBefore);
        }
        if (exportEnabled)
            writeLinearSystemExport(exportPrefix, exportIndex, flag != 0,
                                    initialSolution, exportByTopology);
        return converged;
    }
#endif

    // residual with V=0
    MultPC(b,Z);
    res_o=Dot(Z,b);
    if(res_o==0)
    {
        m_lastRelativeResidual = 0.;
        if (m_collectStats)
        {
            ++m_solveCount;
            std::fprintf(stderr,
                         "XFEMM_PCG_SOLVE solve=%llu iterations=0 "
                         "relative_residual=0 wall_s=%.6f\n",
                         m_solveCount - 1,
                         std::chrono::duration<double>(
                             std::chrono::steady_clock::now() - solveStarted).count());
        }
        if (exportEnabled)
            writeLinearSystemExport(exportPrefix, exportIndex, flag != 0,
                                    initialSolution, exportByTopology);
        return true;
    }

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
        ++m_lastIterations;
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

    m_lastRelativeResidual = er;
    if (m_collectStats)
    {
        const double wallSeconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - solveStarted).count();
        const double spmvSeconds = m_spmvSeconds - spmvSecondsBefore;
        const double preconditionerSeconds =
            m_preconditionerSeconds - preconditionerSecondsBefore;
        const double dotSeconds = m_dotSeconds - dotSecondsBefore;
        const double packSeconds = m_packSeconds - packSecondsBefore;
        const double otherSeconds = std::max(
            0., wallSeconds - packSeconds - spmvSeconds -
                    preconditionerSeconds - dotSeconds);
        ++m_solveCount;
        m_totalIterations += static_cast<unsigned long long>(m_lastIterations);
        std::fprintf(
            stderr,
            "XFEMM_PCG_SOLVE solve=%llu iterations=%d relative_residual=%.9g "
            "wall_s=%.6f pack_s=%.6f spmv_s=%.6f preconditioner_s=%.6f "
            "dot_s=%.6f other_s=%.6f pack_calls=%llu spmv_calls=%llu "
            "preconditioner_calls=%llu dot_calls=%llu\n",
            m_solveCount - 1, m_lastIterations, m_lastRelativeResidual,
            wallSeconds, packSeconds, spmvSeconds, preconditionerSeconds,
            dotSeconds, otherSeconds, m_packCalls - packCallsBefore,
            m_spmvCalls - spmvCallsBefore,
            m_preconditionerCalls - preconditionerCallsBefore,
            m_dotCalls - dotCallsBefore);
    }

    if (exportEnabled)
        writeLinearSystemExport(exportPrefix, exportIndex, flag != 0,
                                initialSolution, exportByTopology);

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
                {
                    ++m_lastIterations;
                    pAp = 0.;
                }
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

    m_lastRelativeResidual = active ? er : 0.;

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
