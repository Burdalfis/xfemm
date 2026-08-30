#ifndef XFEMM_LINSOLVE_CUDSS_LINEAR_SYSTEM_BACKEND_H
#define XFEMM_LINSOLVE_CUDSS_LINEAR_SYSTEM_BACKEND_H

#include "LinearSystemBackend.h"

#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace femm {

/** Runtime probe used by optional GPU regression tests and applications. */
bool cudssDeviceAvailable() noexcept;

/** Immutable AGE superset selected before numeric assembly begins. */
struct CudssBucketDefinition {
    std::string identity;
    /** Additional upper-triangular structural entries retained as numeric zero. */
    std::vector<std::pair<std::int32_t, std::int32_t>> upperEntries;
};

struct CudssBackendOptions {
    bool deterministic = false;
    /** Maximum resident assembly/cuDSS contexts. Must allow committed + trial. */
    std::size_t bucketCacheCapacity = 2;
};

/**
 * Native-real FP64 cuDSS backend.
 *
 * FEM assembly remains on the CPU in CBigLinProb.  Each AGE bucket owns a
 * persistent assembly graph and an immutable cuDSS analysis context.  Values
 * and RHS are updated in place; cuDSS performs numeric refactorization and a
 * direct SPD solve without using solution() as an initial guess.
 */
class CudssLinearSystemBackend final : public LinearSystemBackend<double> {
public:
    explicit CudssLinearSystemBackend(CudssBackendOptions options = {});
    ~CudssLinearSystemBackend() override;

    CudssLinearSystemBackend(const CudssLinearSystemBackend &) = delete;
    CudssLinearSystemBackend &operator=(const CudssLinearSystemBackend &) = delete;

    /** Select a bucket before wipe/assembly. An analyzed bucket is immutable. */
    void activateBucket(const CudssBucketDefinition &definition);
    /** Pin the accepted bucket; updating the pin may evict an older context. */
    void setCommittedBucket(const std::string &identity);
    std::size_t residentBucketCount() const;
    std::size_t bucketEvictionCount() const;

    ScalarType scalar_type() const override { return ScalarType::Real; }
    bool create(int dimension, int bandwidth, int node_count = -1) override;
    int dimension() const override;
    void wipe() override;
    void put(double value, int row, int col, int matrix = 0) override;
    void add_to(double value, int row, int col, int matrix = 0) override;
    void add_symmetric_3x3(std::size_t elementIndex,
                           const int nodes[3],
                           const double values[6]) override;
    double get(int row, int col, int matrix = 0) override;
    void set_value(int i, double x) override;
    void constrain_periodic(int a, int b, bool antiperiodic) override;
    ScalarView<double> &rhs() override;
    const ScalarView<double> &rhs() const override;
    ScalarView<double> &solution() override;
    const ScalarView<double> &solution() const override;
    ScalarView<int> &node_flag() override;
    const ScalarView<int> &node_flag() const override;
    ScalarView<double> &scratch() override;
    bool newton() const override { return false; }
    void set_newton(bool) override {}
    double precision() const override;
    void set_precision(double p) override;
    SolveReport solve(const SolveOptions &options) override;
    RetainedFactorizationSolveReport solveRetainedFactorization(
        const std::vector<double> &rightHandSides,
        std::size_t rightHandSideCount) override;
    void reset_diagnostics() override;
    LinearSystemDiagnostics diagnostics() const override;

private:
    class Implementation;
    std::unique_ptr<Implementation> m_impl;
};

} // namespace femm

#endif
