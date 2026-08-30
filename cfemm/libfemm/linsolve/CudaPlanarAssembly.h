#ifndef XFEMM_LINSOLVE_CUDA_PLANAR_ASSEMBLY_H
#define XFEMM_LINSOLVE_CUDA_PLANAR_ASSEMBLY_H

#include "LinearSystemBackend.h"

#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

namespace femm {

struct CudaPlanarMatrixAdd
{
    std::int32_t row = 0;
    std::int32_t column = 0;
    double value = 0;
};

enum class CudaPlanarConstraintKind { Dirichlet, Periodic, Antiperiodic };

struct CudaPlanarConstraint
{
    CudaPlanarConstraintKind kind = CudaPlanarConstraintKind::Dirichlet;
    std::int32_t a = 0;
    std::int32_t b = 0;
    double value = 0;
};

struct CudaPlanarAssemblyTimings
{
    double clearMs = 0;
    double materialMs = 0;
    double elementMs = 0;
    double scatterMs = 0;
    double ageUploadMs = 0;
    double constraintMs = 0;
    std::uint64_t transferBytes = 0;
};

struct CudaPlanarBhResult
{
    double reluctivity = 0;
    double differentialReluctivity = 0;
};

/** Focused FP64 constitutive-parity hook used by CUDA regression tests. */
std::vector<CudaPlanarBhResult> evaluateCudaPlanarBh(
    const PlanarAssemblyMaterial &material,
    const std::vector<double> &bhFluxDensity,
    const std::vector<double> &bhField,
    const std::vector<double> &bhSlope,
    const std::vector<double> &fluxDensitySamples);

/** Bucket-scoped CUDA numerical assembly plan. The implementation is in a
 * CUDA translation unit so normal host-only builds never see CUDA types. */
class CudaPlanarAssembly
{
public:
    CudaPlanarAssembly(const PlanarAssemblyPlan &plan,
                       const std::vector<std::int32_t> &lowerRows,
                       const std::vector<std::int32_t> &lowerColumns,
                       PlanarAssemblyBackend mode);
    ~CudaPlanarAssembly();
    CudaPlanarAssembly(const CudaPlanarAssembly &) = delete;
    CudaPlanarAssembly &operator=(const CudaPlanarAssembly &) = delete;

    CudaPlanarAssemblyTimings assemble(
        void *stream, double *deviceValues, double *deviceRhs,
        double *deviceSolution, const double *hostSolution,
        const PlanarAssemblyState &state,
        const std::vector<CudaPlanarMatrixAdd> &ageContributions);
    double applyConstraints(void *stream, double *deviceValues,
                            double *deviceRhs,
                            const std::vector<CudaPlanarConstraint> &constraints);
    double relativeResidual(void *stream, const double *deviceValues,
                            const double *deviceRhs,
                            const double *deviceSolution);

    void downloadMatrixAndRhs(void *stream, const double *deviceValues,
                              const double *deviceRhs,
                              std::vector<double> &values,
                              std::vector<double> &rhs) const;
    /** Storage owned by the bucket-scoped assembly plan, excluding the CSR,
     * RHS, solution, and cuDSS-internal allocations owned by CudssContext. */
    std::uint64_t deviceBytes() const;

private:
    class Implementation;
    std::unique_ptr<Implementation> m_impl;
};

} // namespace femm

#endif
