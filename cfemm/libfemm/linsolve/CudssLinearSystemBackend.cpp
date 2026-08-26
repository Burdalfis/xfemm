#include "CudssLinearSystemBackend.h"

#include "spars.h"

#include <cudss.h>
#include <cuda_runtime.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>
#include <stdexcept>
#include <unordered_map>

namespace femm {
bool cudssDeviceAvailable() noexcept
{
    int count = 0;
    return cudaGetDeviceCount(&count) == cudaSuccess && count > 0;
}

namespace {

using Clock = std::chrono::steady_clock;

double elapsedMs(Clock::time_point start)
{
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

void cudaCheck(cudaError_t status, const char *operation)
{
    if (status != cudaSuccess)
        throw std::runtime_error(std::string(operation) + ": " +
                                 cudaGetErrorString(status));
}

void cudssCheck(cudssStatus_t status, const char *operation)
{
    if (status != CUDSS_STATUS_SUCCESS)
        throw std::runtime_error(std::string(operation) + ": cuDSS status " +
                                 std::to_string(static_cast<int>(status)));
}

#define XFEMM_CUDA_CHECK(expr) cudaCheck((expr), #expr)
#define XFEMM_CUDSS_CHECK(expr) cudssCheck((expr), #expr)

template <typename Function>
double timeStream(cudaStream_t stream, Function &&function)
{
    cudaEvent_t begin = nullptr;
    cudaEvent_t end = nullptr;
    XFEMM_CUDA_CHECK(cudaEventCreate(&begin));
    XFEMM_CUDA_CHECK(cudaEventCreate(&end));
    try {
        XFEMM_CUDA_CHECK(cudaEventRecord(begin, stream));
        function();
        XFEMM_CUDA_CHECK(cudaEventRecord(end, stream));
        XFEMM_CUDA_CHECK(cudaEventSynchronize(end));
        float milliseconds = 0;
        XFEMM_CUDA_CHECK(cudaEventElapsedTime(&milliseconds, begin, end));
        cudaEventDestroy(end);
        cudaEventDestroy(begin);
        return milliseconds;
    } catch (...) {
        if (end) cudaEventDestroy(end);
        if (begin) cudaEventDestroy(begin);
        throw;
    }
}

template <typename T>
class DeviceBuffer {
public:
    DeviceBuffer() = default;
    explicit DeviceBuffer(std::size_t count) { allocate(count); }
    DeviceBuffer(const DeviceBuffer &) = delete;
    DeviceBuffer &operator=(const DeviceBuffer &) = delete;
    ~DeviceBuffer() { if (m_data) cudaFree(m_data); }

    void allocate(std::size_t count)
    {
        if (m_data) throw std::logic_error("CUDA buffer is already allocated");
        m_count = count;
        if (count) XFEMM_CUDA_CHECK(cudaMalloc(&m_data, count * sizeof(T)));
    }
    T *get() const { return m_data; }
    std::size_t bytes() const { return m_count * sizeof(T); }

private:
    T *m_data = nullptr;
    std::size_t m_count = 0;
};

std::uint64_t topologyHash(const std::vector<std::int32_t> &rows,
                           const std::vector<std::int32_t> &columns)
{
    std::uint64_t hash = 1469598103934665603ull;
    const auto mix = [&hash](const void *pointer, std::size_t size) {
        const auto *bytes = static_cast<const unsigned char *>(pointer);
        for (std::size_t i = 0; i < size; ++i) {
            hash ^= bytes[i];
            hash *= 1099511628211ull;
        }
    };
    mix(rows.data(), rows.size() * sizeof(rows[0]));
    mix(columns.data(), columns.size() * sizeof(columns[0]));
    return hash;
}

struct HostCsr {
    std::vector<std::int32_t> rows;
    std::vector<std::int32_t> columns;
    std::vector<double> values;
};

HostCsr transposeUpperToLower(const HostCsr &upper)
{
    if (upper.rows.empty()) return {};
    const std::size_t dimension = upper.rows.size() - 1;
    HostCsr lower;
    lower.rows.assign(dimension + 1, 0);
    for (std::int32_t column : upper.columns) {
        if (column < 0 || static_cast<std::size_t>(column) >= dimension)
            throw std::out_of_range("upper CSR column is outside matrix bounds");
        ++lower.rows[static_cast<std::size_t>(column) + 1];
    }
    for (std::size_t row = 0; row < dimension; ++row)
        lower.rows[row + 1] += lower.rows[row];
    lower.columns.resize(upper.columns.size());
    lower.values.resize(upper.values.size());
    std::vector<std::int32_t> next = lower.rows;
    for (std::size_t row = 0; row < dimension; ++row) {
        for (std::int32_t j = upper.rows[row]; j < upper.rows[row + 1]; ++j) {
            const std::size_t destinationRow =
                static_cast<std::size_t>(upper.columns[static_cast<std::size_t>(j)]);
            const std::int32_t destination = next[destinationRow]++;
            lower.columns[static_cast<std::size_t>(destination)] =
                static_cast<std::int32_t>(row);
            lower.values[static_cast<std::size_t>(destination)] =
                upper.values[static_cast<std::size_t>(j)];
        }
    }
    return lower;
}

HostCsr makeUnion(const HostCsr &exact,
                  const std::vector<std::pair<std::int32_t, std::int32_t>> &extra,
                  int dimension)
{
    std::vector<std::vector<std::int32_t>> additions(static_cast<std::size_t>(dimension));
    for (auto entry : extra) {
        if (entry.second < entry.first) std::swap(entry.first, entry.second);
        if (entry.first < 0 || entry.second >= dimension)
            throw std::out_of_range("cuDSS bucket entry is outside matrix bounds");
        additions[static_cast<std::size_t>(entry.first)].push_back(entry.second);
    }
    for (auto &row : additions) {
        std::sort(row.begin(), row.end());
        row.erase(std::unique(row.begin(), row.end()), row.end());
    }

    HostCsr result;
    result.rows.resize(static_cast<std::size_t>(dimension) + 1);
    for (int row = 0; row < dimension; ++row) {
        result.rows[static_cast<std::size_t>(row)] =
            static_cast<std::int32_t>(result.columns.size());
        const std::int32_t begin = exact.rows[static_cast<std::size_t>(row)];
        const std::int32_t end = exact.rows[static_cast<std::size_t>(row) + 1];
        const auto &added = additions[static_cast<std::size_t>(row)];
        auto add = added.begin();
        std::int32_t current = begin;
        while (current < end || add != added.end()) {
            const std::int32_t exactColumn = current < end
                ? exact.columns[static_cast<std::size_t>(current)]
                : std::numeric_limits<std::int32_t>::max();
            const std::int32_t addedColumn = add != added.end()
                ? *add : std::numeric_limits<std::int32_t>::max();
            const std::int32_t column = std::min(exactColumn, addedColumn);
            result.columns.push_back(column);
            result.values.push_back(exactColumn == column
                ? exact.values[static_cast<std::size_t>(current)] : 0.0);
            if (exactColumn == column) ++current;
            if (addedColumn == column) ++add;
        }
    }
    result.rows.back() = static_cast<std::int32_t>(result.columns.size());
    return result;
}

bool scatterValues(const HostCsr &exact, const std::vector<std::int32_t> &targetRows,
                   const std::vector<std::int32_t> &targetColumns,
                   std::vector<double> &targetValues)
{
    std::fill(targetValues.begin(), targetValues.end(), 0.0);
    const std::size_t dimension = exact.rows.size() - 1;
    for (std::size_t row = 0; row < dimension; ++row) {
        std::int32_t target = targetRows[row];
        const std::int32_t targetEnd = targetRows[row + 1];
        for (std::int32_t source = exact.rows[row]; source < exact.rows[row + 1]; ++source) {
            const std::int32_t column = exact.columns[static_cast<std::size_t>(source)];
            while (target < targetEnd &&
                   targetColumns[static_cast<std::size_t>(target)] < column)
                ++target;
            if (target == targetEnd ||
                targetColumns[static_cast<std::size_t>(target)] != column)
                return false;
            targetValues[static_cast<std::size_t>(target)] =
                exact.values[static_cast<std::size_t>(source)];
        }
    }
    return true;
}

double relativeResidual(const HostCsr &matrix, const double *rhs,
                        const double *solution)
{
    long double residualSquared = 0;
    long double rhsSquared = 0;
    std::vector<long double> product(matrix.rows.size() - 1, 0.0L);
    for (std::size_t row = 0; row + 1 < matrix.rows.size(); ++row) {
        for (std::int32_t j = matrix.rows[row]; j < matrix.rows[row + 1]; ++j) {
            const std::size_t column =
                static_cast<std::size_t>(matrix.columns[static_cast<std::size_t>(j)]);
            const long double value = matrix.values[static_cast<std::size_t>(j)];
            product[row] += value * solution[column];
            if (column != row) product[column] += value * solution[row];
        }
    }
    for (std::size_t i = 0; i < product.size(); ++i) {
        const long double residual = product[i] - rhs[i];
        residualSquared += residual * residual;
        rhsSquared += static_cast<long double>(rhs[i]) * rhs[i];
    }
    if (rhsSquared == 0) return std::sqrt(static_cast<double>(residualSquared));
    return std::sqrt(static_cast<double>(residualSquared / rhsSquared));
}

class CudssContext {
public:
    CudssContext(const HostCsr &matrix, const double *rhs, bool deterministic,
                 LinearSystemDiagnostics &diagnostics)
        : m_rows(matrix.rows), m_columns(matrix.columns), m_values(matrix.values)
    {
        const auto constructionStarted = Clock::now();
        m_dRows.allocate(m_rows.size());
        m_dColumns.allocate(m_columns.size());
        m_dValues.allocate(m_values.size());
        m_dRhs.allocate(m_rows.size() - 1);
        m_dSolution.allocate(m_rows.size() - 1);
        XFEMM_CUDSS_CHECK(cudssCreate(&m_handle));
        XFEMM_CUDA_CHECK(cudaStreamCreateWithFlags(&m_stream, cudaStreamNonBlocking));
        XFEMM_CUDSS_CHECK(cudssSetStream(m_handle, m_stream));
        XFEMM_CUDSS_CHECK(cudssConfigCreate(&m_config));
        if (deterministic) {
            const int enabled = 1;
            XFEMM_CUDSS_CHECK(cudssConfigSet(m_config, CUDSS_CONFIG_DETERMINISTIC_MODE,
                                             &enabled, sizeof(enabled)));
        }
        XFEMM_CUDSS_CHECK(cudssDataCreate(m_handle, &m_data));
        const std::int64_t n = static_cast<std::int64_t>(m_rows.size() - 1);
        const std::int64_t nnz = static_cast<std::int64_t>(m_columns.size());
        XFEMM_CUDSS_CHECK(cudssMatrixCreateCsr(
            &m_matrix, n, n, nnz, m_dRows.get(), nullptr, m_dColumns.get(),
            m_dValues.get(), CUDSS_R_32I, CUDSS_R_32I, CUDSS_R_64F,
            CUDSS_MTYPE_SPD, CUDSS_MVIEW_LOWER, CUDSS_BASE_ZERO));
        XFEMM_CUDSS_CHECK(cudssMatrixCreateDn(&m_rhs, n, 1, n, m_dRhs.get(),
                                               CUDSS_R_64F, CUDSS_LAYOUT_COL_MAJOR));
        XFEMM_CUDSS_CHECK(cudssMatrixCreateDn(&m_solution, n, 1, n,
                                               m_dSolution.get(), CUDSS_R_64F,
                                               CUDSS_LAYOUT_COL_MAJOR));
        diagnostics.bucketConstructionMs += elapsedMs(constructionStarted);
        diagnostics.hostToDeviceMs += upload(rhs, diagnostics);
        const double analysis = timeStream(m_stream, [&] {
            XFEMM_CUDSS_CHECK(cudssExecute(m_handle, CUDSS_PHASE_ANALYSIS, m_config,
                                           m_data, m_matrix, m_solution, m_rhs));
        });
        diagnostics.symbolicAnalysisMs += analysis;
        std::size_t written = 0;
        if (cudssDataGet(m_handle, m_data, CUDSS_DATA_MEMORY_ESTIMATES,
                         m_memory.data(), m_memory.size() * sizeof(m_memory[0]),
                         &written) != CUDSS_STATUS_SUCCESS)
            m_memory.fill(0);
    }

    ~CudssContext()
    {
        if (m_solution) cudssMatrixDestroy(m_solution);
        if (m_rhs) cudssMatrixDestroy(m_rhs);
        if (m_matrix) cudssMatrixDestroy(m_matrix);
        if (m_data) cudssDataDestroy(m_handle, m_data);
        if (m_config) cudssConfigDestroy(m_config);
        if (m_handle) cudssDestroy(m_handle);
        if (m_stream) cudaStreamDestroy(m_stream);
    }

    bool contains(const HostCsr &exact, LinearSystemDiagnostics &diagnostics) const
    {
        const auto started = Clock::now();
        std::vector<double> ignored(m_values.size());
        const bool result = scatterValues(exact, m_rows, m_columns, ignored);
        diagnostics.sparsePackingMs += elapsedMs(started);
        return result;
    }

    void solve(const HostCsr &exact, const double *rhs, double *solution,
               LinearSystemDiagnostics &diagnostics)
    {
        const auto scatterStarted = Clock::now();
        if (!scatterValues(exact, m_rows, m_columns, m_values))
            throw std::logic_error("attempted to mutate an analyzed cuDSS sparsity graph");
        diagnostics.sparsePackingMs += elapsedMs(scatterStarted);
        diagnostics.hostToDeviceMs += upload(rhs, diagnostics);
        const int phase = m_factorized ? CUDSS_PHASE_REFACTORIZATION
                                       : CUDSS_PHASE_FACTORIZATION;
        diagnostics.numericFactorizationMs += timeStream(m_stream, [&] {
            XFEMM_CUDSS_CHECK(cudssExecute(m_handle, phase, m_config, m_data,
                                           m_matrix, m_solution, m_rhs));
        });
        diagnostics.solveMs += timeStream(m_stream, [&] {
            XFEMM_CUDSS_CHECK(cudssExecute(m_handle, CUDSS_PHASE_SOLVE, m_config,
                                           m_data, m_matrix, m_solution, m_rhs));
        });
        diagnostics.deviceToHostMs += timeStream(m_stream, [&] {
            XFEMM_CUDA_CHECK(cudaMemcpyAsync(
                solution, m_dSolution.get(), m_dSolution.bytes(),
                cudaMemcpyDeviceToHost, m_stream));
        });
        m_factorized = true;

        std::int64_t factorNnz = 0;
        std::size_t written = 0;
        if (cudssDataGet(m_handle, m_data, CUDSS_DATA_LU_NNZ, &factorNnz,
                         sizeof(factorNnz), &written) == CUDSS_STATUS_SUCCESS)
            diagnostics.factorNonzeros =
                static_cast<std::uint64_t>(std::max<std::int64_t>(0, factorNnz));
        diagnostics.matrixNonzeros = m_columns.size();
        diagnostics.permanentDeviceBytes =
            static_cast<std::uint64_t>(std::max<std::int64_t>(0, m_memory[0]));
        diagnostics.peakDeviceBytes =
            static_cast<std::uint64_t>(std::max<std::int64_t>(0, m_memory[1]));
        diagnostics.factorizationRetained = true;
    }

    bool factorized() const { return m_factorized; }

private:
    double upload(const double *rhs, LinearSystemDiagnostics &diagnostics)
    {
        const std::uint64_t bytes = m_dValues.bytes() + m_dRhs.bytes();
        diagnostics.hostToDeviceBytes += bytes;
        return timeStream(m_stream, [&] {
            if (!m_structureUploaded) {
                XFEMM_CUDA_CHECK(cudaMemcpyAsync(m_dRows.get(), m_rows.data(),
                                                  m_dRows.bytes(), cudaMemcpyHostToDevice,
                                                  m_stream));
                XFEMM_CUDA_CHECK(cudaMemcpyAsync(m_dColumns.get(), m_columns.data(),
                                                  m_dColumns.bytes(), cudaMemcpyHostToDevice,
                                                  m_stream));
                diagnostics.hostToDeviceBytes += m_dRows.bytes() + m_dColumns.bytes();
                XFEMM_CUDA_CHECK(cudaMemsetAsync(m_dSolution.get(), 0,
                                                  m_dSolution.bytes(), m_stream));
                m_structureUploaded = true;
            }
            XFEMM_CUDA_CHECK(cudaMemcpyAsync(m_dValues.get(), m_values.data(),
                                              m_dValues.bytes(), cudaMemcpyHostToDevice,
                                              m_stream));
            XFEMM_CUDA_CHECK(cudaMemcpyAsync(m_dRhs.get(), rhs, m_dRhs.bytes(),
                                              cudaMemcpyHostToDevice, m_stream));
        });
    }

    std::vector<std::int32_t> m_rows;
    std::vector<std::int32_t> m_columns;
    std::vector<double> m_values;
    DeviceBuffer<std::int32_t> m_dRows;
    DeviceBuffer<std::int32_t> m_dColumns;
    DeviceBuffer<double> m_dValues;
    DeviceBuffer<double> m_dRhs;
    DeviceBuffer<double> m_dSolution;
    cudssHandle_t m_handle = nullptr;
    cudaStream_t m_stream = nullptr;
    cudssConfig_t m_config = nullptr;
    cudssData_t m_data = nullptr;
    cudssMatrix_t m_matrix = nullptr;
    cudssMatrix_t m_rhs = nullptr;
    cudssMatrix_t m_solution = nullptr;
    std::array<std::int64_t, 16> m_memory{};
    bool m_structureUploaded = false;
    bool m_factorized = false;
};

} // namespace

class CudssLinearSystemBackend::Implementation {
public:
    struct Bucket {
        std::unique_ptr<CBigLinProb> assembly;
        CudssBucketDefinition definition;
        std::unique_ptr<CudssContext> unionContext;
        std::unordered_map<std::uint64_t, std::unique_ptr<CudssContext>> exactContexts;
    };

    explicit Implementation(CudssBackendOptions requested) : options(requested)
    {
        diagnostics.deterministic = options.deterministic;
    }

    Bucket &ensureActive()
    {
        if (active) return *active;
        activate({"exact-default", {}});
        return *active;
    }

    void bind(Bucket &bucket)
    {
        rhs.assign(bucket.assembly->b, static_cast<std::size_t>(dimension));
        solution.assign(bucket.assembly->V, static_cast<std::size_t>(dimension));
        scratch.assign(bucket.assembly->P, static_cast<std::size_t>(dimension));
        flags.assign(bucket.assembly->Q, static_cast<std::size_t>(dimension));
    }

    void activate(CudssBucketDefinition definition)
    {
        const auto lookupStarted = Clock::now();
        if (dimension <= 0)
            throw std::logic_error("create must be called before selecting a cuDSS bucket");
        if (definition.identity.empty()) definition.identity = "exact-default";
        if (active && active->definition.identity == definition.identity) {
            diagnostics.bucketIdentity = definition.identity;
            diagnostics.bucketReused = true;
            diagnostics.bucketLookupMs += elapsedMs(lookupStarted);
            return;
        }

        auto found = buckets.find(definition.identity);
        diagnostics.bucketLookupMs += elapsedMs(lookupStarted);
        if (found == buckets.end()) {
            const auto constructionStarted = Clock::now();
            auto bucket = std::make_unique<Bucket>();
            std::sort(definition.upperEntries.begin(), definition.upperEntries.end());
            definition.upperEntries.erase(
                std::unique(definition.upperEntries.begin(), definition.upperEntries.end()),
                definition.upperEntries.end());
            bucket->definition = std::move(definition);
            bucket->assembly = std::make_unique<CBigLinProb>();
            if (!bucket->assembly->Create(dimension, bandwidth))
                throw std::runtime_error("could not allocate cuDSS bucket assembly matrix");
            // Seed the retained CPU assembly graph with the same immutable AGE
            // union that will be analyzed by cuDSS. CBigLinProb::Wipe preserves
            // structural entries, so subsequent angle-specific assembly only
            // changes numeric values and cannot grow the graph for an expected
            // topology in this bucket.
            for (const auto &entry : bucket->definition.upperEntries)
                bucket->assembly->Put(0.0, entry.first, entry.second);
            found = buckets.emplace(bucket->definition.identity, std::move(bucket)).first;
            diagnostics.bucketConstructionMs += elapsedMs(constructionStarted);
            diagnostics.bucketReused = false;
        } else {
            diagnostics.bucketReused = true;
        }
        const auto switchStarted = Clock::now();
        if (active)
            sessionSolution.assign(active->assembly->V,
                                   active->assembly->V + dimension);
        active = found->second.get();
        if (!sessionSolution.empty())
            std::copy(sessionSolution.begin(), sessionSolution.end(), active->assembly->V);
        bind(*active);
        diagnostics.bucketIdentity = active->definition.identity;
        diagnostics.bucketSwitchMs += elapsedMs(switchStarted);
    }

    HostCsr exactMatrix()
    {
        HostCsr exact;
        active->assembly->copyUpperCsr(exact.rows, exact.columns, exact.values);
        return exact;
    }

    CudssBackendOptions options;
    int dimension = 0;
    int bandwidth = 0;
    double precision = 1e-8;
    std::map<std::string, std::unique_ptr<Bucket>> buckets;
    Bucket *active = nullptr;
    std::vector<double> sessionSolution;
    ScalarView<double> rhs;
    ScalarView<double> solution;
    ScalarView<double> scratch;
    ScalarView<int> flags;
    LinearSystemDiagnostics diagnostics;
};

CudssLinearSystemBackend::CudssLinearSystemBackend(CudssBackendOptions options)
    : m_impl(std::make_unique<Implementation>(options)) {}

CudssLinearSystemBackend::~CudssLinearSystemBackend() = default;

void CudssLinearSystemBackend::activateBucket(const CudssBucketDefinition &definition)
{
    m_impl->activate(definition);
}

bool CudssLinearSystemBackend::create(int dimension, int bandwidth, int)
{
    if (dimension <= 0 || m_impl->dimension != 0) return false;
    m_impl->dimension = dimension;
    m_impl->bandwidth = bandwidth;
    m_impl->sessionSolution.assign(static_cast<std::size_t>(dimension), 0.0);
    return true;
}

int CudssLinearSystemBackend::dimension() const { return m_impl->dimension; }
void CudssLinearSystemBackend::wipe() { m_impl->ensureActive().assembly->Wipe(); }
void CudssLinearSystemBackend::put(double v, int r, int c, int)
{ m_impl->ensureActive().assembly->Put(v, r, c); }
void CudssLinearSystemBackend::add_to(double v, int r, int c, int)
{ m_impl->ensureActive().assembly->AddTo(v, r, c); }
double CudssLinearSystemBackend::get(int r, int c, int)
{ return m_impl->ensureActive().assembly->Get(r, c); }
void CudssLinearSystemBackend::set_value(int i, double x)
{ m_impl->ensureActive().assembly->SetValue(i, x); }
void CudssLinearSystemBackend::constrain_periodic(int a, int b, bool anti)
{
    if (anti) m_impl->ensureActive().assembly->AntiPeriodicity(a, b);
    else m_impl->ensureActive().assembly->Periodicity(a, b);
}
ScalarView<double> &CudssLinearSystemBackend::rhs()
{ m_impl->ensureActive(); return m_impl->rhs; }
const ScalarView<double> &CudssLinearSystemBackend::rhs() const
{ return const_cast<CudssLinearSystemBackend *>(this)->rhs(); }
ScalarView<double> &CudssLinearSystemBackend::solution()
{ m_impl->ensureActive(); return m_impl->solution; }
const ScalarView<double> &CudssLinearSystemBackend::solution() const
{ return const_cast<CudssLinearSystemBackend *>(this)->solution(); }
ScalarView<int> &CudssLinearSystemBackend::node_flag()
{ m_impl->ensureActive(); return m_impl->flags; }
const ScalarView<int> &CudssLinearSystemBackend::node_flag() const
{ return const_cast<CudssLinearSystemBackend *>(this)->node_flag(); }
ScalarView<double> &CudssLinearSystemBackend::scratch()
{ m_impl->ensureActive(); return m_impl->scratch; }
double CudssLinearSystemBackend::precision() const { return m_impl->precision; }
void CudssLinearSystemBackend::set_precision(double p)
{
    m_impl->precision = p;
    for (auto &entry : m_impl->buckets) entry.second->assembly->Precision = p;
}

SolveReport CudssLinearSystemBackend::solve(const SolveOptions &options)
{
    auto &bucket = m_impl->ensureActive();
    bucket.assembly->Precision = options.tolerance > 0 ? options.tolerance
                                                       : m_impl->precision;
    const auto packingStarted = Clock::now();
    HostCsr exactUpper = m_impl->exactMatrix();
    HostCsr exact = transposeUpperToLower(exactUpper);
    m_impl->diagnostics.sparsePackingMs += elapsedMs(packingStarted);
    CudssContext *context = nullptr;
    if (!bucket.unionContext) {
        HostCsr united = transposeUpperToLower(
            makeUnion(exactUpper, bucket.definition.upperEntries,
                      m_impl->dimension));
        bucket.unionContext = std::make_unique<CudssContext>(
            united, bucket.assembly->b, m_impl->options.deterministic,
            m_impl->diagnostics);
        context = bucket.unionContext.get();
        m_impl->diagnostics.symbolicReused = false;
    } else if (bucket.unionContext->contains(exact, m_impl->diagnostics)) {
        context = bucket.unionContext.get();
        m_impl->diagnostics.symbolicReused = true;
    } else {
        // Never widen an analyzed graph. Cache the unforeseen exact graph as
        // an explicitly diagnosed safe fallback context.
        m_impl->diagnostics.exactTopologyFallback = true;
        const std::uint64_t hash = topologyHash(exact.rows, exact.columns);
        auto found = bucket.exactContexts.find(hash);
        if (found == bucket.exactContexts.end()) {
            auto inserted = bucket.exactContexts.emplace(
                hash, std::make_unique<CudssContext>(
                    exact, bucket.assembly->b, m_impl->options.deterministic,
                    m_impl->diagnostics));
            found = inserted.first;
            m_impl->diagnostics.symbolicReused = false;
        } else {
            m_impl->diagnostics.symbolicReused = true;
        }
        context = found->second.get();
    }

    context->solve(exact, bucket.assembly->b, bucket.assembly->V,
                   m_impl->diagnostics);
    ++m_impl->diagnostics.linearSolves;
    m_impl->sessionSolution.assign(bucket.assembly->V,
                                   bucket.assembly->V + m_impl->dimension);

    SolveReport report;
    report.solver = "cudss-direct-spd-fp64";
    report.iterations = 1;
    const auto residualStarted = Clock::now();
    report.relative_residual = relativeResidual(
        exact, bucket.assembly->b, bucket.assembly->V);
    m_impl->diagnostics.residualEvaluationMs += elapsedMs(residualStarted);
    report.converged = std::isfinite(report.relative_residual) &&
                       report.relative_residual <= bucket.assembly->Precision;
    m_impl->diagnostics.solver = report.solver;
    m_impl->diagnostics.lastRelativeResidual = report.relative_residual;
    m_impl->diagnostics.allConverged =
        m_impl->diagnostics.allConverged && report.converged;
    if (!report.converged)
        throw std::runtime_error(
            "cuDSS conventional residual " +
            std::to_string(report.relative_residual) +
            " exceeds tolerance " + std::to_string(bucket.assembly->Precision) +
            " for dimension " + std::to_string(m_impl->dimension) +
            " and stored triangular CSR NNZ " +
            std::to_string(exact.columns.size()));
    return report;
}

void CudssLinearSystemBackend::reset_diagnostics()
{
    const bool deterministic = m_impl->options.deterministic;
    const std::string bucket = m_impl->active
        ? m_impl->active->definition.identity : std::string();
    m_impl->diagnostics = {};
    m_impl->diagnostics.deterministic = deterministic;
    m_impl->diagnostics.bucketIdentity = bucket;
    m_impl->diagnostics.allConverged = true;
}

LinearSystemDiagnostics CudssLinearSystemBackend::diagnostics() const
{ return m_impl->diagnostics; }

} // namespace femm
