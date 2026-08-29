#include "PersistentMotorSession.h"

#include "FemmReader.h"
#include "linsolve/CudssLinearSystemBackend.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

std::shared_ptr<femm::FemmProblem> loadProblem(const std::string &path)
{
    auto problem = std::make_shared<femm::FemmProblem>(
        femm::FileType::MagneticsFile);
    std::ostringstream errors;
    {
        femm::MagneticsReader reader(problem, errors);
        if (reader.parse(path) != femm::F_FILE_OK)
            throw std::runtime_error("could not parse motor fixture: " + errors.str());
    }
    // Use a tighter nonlinear/Legacy-PCG tolerance for the backend parity
    // oracle. At the fixture's production 1e-8 tolerance, Legacy truncation
    // dominates tiny flux and zero-order AGE quantities even though cuDSS's
    // conventional residual is near machine precision.
    problem->Precision = 1e-10;
    return problem;
}

double relativeVectorError(const std::vector<double> &a,
                           const std::vector<double> &b)
{
    if (a.size() != b.size()) return std::numeric_limits<double>::infinity();
    long double error = 0;
    long double reference = 0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        const long double delta = static_cast<long double>(a[i]) - b[i];
        error += delta * delta;
        reference += static_cast<long double>(a[i]) * a[i];
    }
    return std::sqrt(static_cast<double>(error /
        std::max<long double>(reference, std::numeric_limits<double>::min())));
}

double maximumVectorError(const std::vector<double> &a,
                          const std::vector<double> &b)
{
    if (a.size() != b.size()) return std::numeric_limits<double>::infinity();
    double result = 0;
    for (std::size_t i = 0; i < a.size(); ++i)
        result = std::max(result, std::abs(a[i] - b[i]));
    return result;
}

bool close(double a, double b, double relative, double absolute)
{
    return std::abs(a - b) <= absolute +
           relative * std::max(std::abs(a), std::abs(b));
}

bool close(const CComplex &a, const CComplex &b,
           double relative, double absolute)
{
    return close(a.re, b.re, relative, absolute) &&
           close(a.im, b.im, relative, absolute);
}

int fail(const std::string &message)
{
    std::cerr << message << '\n';
    return 1;
}

} // namespace

int main(int argc, char **argv)
{
    if (argc != 2) return fail("expected a sliding-AGE motor .fem path");
    if (!femm::cudssDeviceAvailable()) {
        std::cerr << "No CUDA device; skipping persistent cuDSS regression\n";
        return 77;
    }

    try {
        femm::PersistentMotorSession cpu(
            femm::ModelDefinition(loadProblem(argv[1])));
        femm::PersistentMotorSession gpu(
            femm::ModelDefinition(loadProblem(argv[1])),
            femm::CudssSessionOptions{1.2, false, false, 2});

        gpu.initialize();
        const auto initialization = gpu.initializationDiagnostics();
        if (initialization.sessionInitializationMs <= 0 ||
            initialization.bucketConstructionMs <= 0 ||
            initialization.symbolicAnalysisMs <= 0 ||
            !initialization.factorizationRetained ||
            initialization.exactTopologyFallback)
            return fail("cuDSS session initialization did not prime the initial bucket");
        if (initialization.matrixNonzeros == 0 ||
            initialization.factorNonzeros == 0 ||
            initialization.serializationPostprocessorMs <= 0)
            return fail("cuDSS initialization/resource diagnostics were not populated");

        const auto cpuReference = cpu.evaluateTrial();
        const auto gpuTrial = gpu.evaluateTrial();
        if (!gpuTrial.real || !cpuReference.real)
            return fail("real solution payload missing");
        if (!gpuTrial.diagnostics.converged ||
            !std::isfinite(gpuTrial.diagnostics.relativeResidual) ||
            gpuTrial.diagnostics.relativeResidual >
                gpu.analysis().model().problem().Precision)
            return fail("cuDSS residual/convergence gate failed");
        if (gpuTrial.diagnostics.symbolicAnalysisMs != 0 ||
            gpuTrial.diagnostics.bucketConstructionMs != 0 ||
            !gpuTrial.diagnostics.symbolicReused ||
            !gpuTrial.diagnostics.bucketReused ||
            !gpuTrial.diagnostics.factorizationRetained)
            return fail("hot evaluation did not reuse initialized cuDSS resources");
        if (gpuTrial.diagnostics.exactTopologyFallback)
            return fail("measured narrow AGE bucket unexpectedly missed its initial graph");
        if (gpuTrial.diagnostics.matrixNonzeros == 0 ||
            gpuTrial.diagnostics.factorNonzeros == 0)
            return fail("cuDSS sparse diagnostics were not populated");
        if (!gpuTrial.diagnostics.nonlinearWarmStartUsed ||
            gpuTrial.diagnostics.nonlinearIterations <= 0 ||
            cpuReference.diagnostics.nonlinearIterations <= 0 ||
            gpuTrial.diagnostics.nonlinearIterations >
                cpuReference.diagnostics.nonlinearIterations)
            return fail("nonlinear initialization/warm-start convergence regression");
        if (gpuTrial.diagnostics.linearSolver != "cudss-direct-spd-fp64" ||
            gpuTrial.diagnostics.legacyFallback)
            return fail("cuDSS evaluation unexpectedly used another linear backend");

        const double nodalError = relativeVectorError(
            cpuReference.real->nodal.magneticVectorPotential,
            gpuTrial.real->nodal.magneticVectorPotential);
        const double nodalMaxError = maximumVectorError(
            cpuReference.real->nodal.magneticVectorPotential,
            gpuTrial.real->nodal.magneticVectorPotential);
        std::cout << std::setprecision(17)
                  << "persistent_session_initial"
                  << " cpu_torque=" << cpuReference.real->torque
                  << " gpu_torque=" << gpuTrial.real->torque
                  << " torque_abs_error="
                  << std::abs(cpuReference.real->torque - gpuTrial.real->torque)
                  << " cpu_nodes="
                  << cpuReference.real->nodal.magneticVectorPotential.size()
                  << " gpu_nodes="
                  << gpuTrial.real->nodal.magneticVectorPotential.size()
                  << " nodal_relative_error=" << nodalError
                  << " nodal_max_error=" << nodalMaxError
                  << " gpu_relative_residual="
                  << gpuTrial.diagnostics.relativeResidual
                  << " cpu_newton_iterations="
                  << cpuReference.diagnostics.nonlinearIterations
                  << " gpu_newton_iterations="
                  << gpuTrial.diagnostics.nonlinearIterations << '\n';
        if (nodalError > 2e-6 || nodalMaxError > 1e-8)
            return fail("CPU/cuDSS nodal solution regression: relative=" +
                        std::to_string(nodalError) + " max=" +
                        std::to_string(nodalMaxError));
        if (!close(cpuReference.real->torque, gpuTrial.real->torque, 5e-6, 1e-10))
            return fail("CPU/cuDSS torque regression: cpu=" +
                        std::to_string(cpuReference.real->torque) + " gpu=" +
                        std::to_string(gpuTrial.real->torque));
        if (gpuTrial.real->circuits.size() != 3 ||
            cpuReference.real->circuits.size() != 3)
            return fail("three circuit flux linkages were not returned");
        for (std::size_t i = 0; i < 3; ++i) {
            std::cout << std::setprecision(17)
                      << "persistent_session_flux circuit=" << i
                      << " cpu=" << cpuReference.real->circuits[i].fluxLinkage
                      << " gpu=" << gpuTrial.real->circuits[i].fluxLinkage
                      << " abs_error="
                      << std::abs(cpuReference.real->circuits[i].fluxLinkage -
                                  gpuTrial.real->circuits[i].fluxLinkage) << '\n';
            if (!close(cpuReference.real->circuits[i].fluxLinkage,
                       gpuTrial.real->circuits[i].fluxLinkage, 5e-5, 5e-11))
                return fail("CPU/cuDSS circuit linkage regression");
        }
        if (cpuReference.real->airGaps.size() != gpuTrial.real->airGaps.size() ||
            gpuTrial.real->airGaps.empty())
            return fail("AGE result payload mismatch");
        for (std::size_t g = 0; g < gpuTrial.real->airGaps.size(); ++g) {
            const auto &a = cpuReference.real->airGaps[g];
            const auto &b = gpuTrial.real->airGaps[g];
            if (a.name != b.name || a.harmonics.size() != b.harmonics.size() ||
                !close(a.torque, b.torque, 5e-6, 1e-10) ||
                !close(a.centerVectorPotential, b.centerVectorPotential,
                       5e-6, 5e-9)) {
                std::ostringstream message;
                message << std::setprecision(17)
                        << "AGE quantity/count regression gap=" << g
                        << " cpu_name=" << a.name << " gpu_name=" << b.name
                        << " cpu_count=" << a.harmonics.size()
                        << " gpu_count=" << b.harmonics.size()
                        << " cpu_torque=" << a.torque
                        << " gpu_torque=" << b.torque
                        << " cpu_aco=(" << a.centerVectorPotential.re << ','
                        << a.centerVectorPotential.im << ')'
                        << " gpu_aco=(" << b.centerVectorPotential.re << ','
                        << b.centerVectorPotential.im << ')';
                return fail(message.str());
            }
            for (std::size_t h = 0; h < a.harmonics.size(); ++h) {
                const auto &ah = a.harmonics[h];
                const auto &bh = b.harmonics[h];
                if (ah.order != bh.order ||
                    !close(ah.radialCos, bh.radialCos, 5e-6, 1e-9) ||
                    !close(ah.radialSin, bh.radialSin, 5e-6, 1e-9) ||
                    !close(ah.tangentialCos, bh.tangentialCos, 5e-6, 1e-9) ||
                    !close(ah.tangentialSin, bh.tangentialSin, 5e-6, 1e-9)) {
                    std::ostringstream message;
                    message << std::setprecision(17)
                            << "AGE harmonic regression gap=" << g
                            << " harmonic=" << h
                            << " order_cpu=" << ah.order
                            << " order_gpu=" << bh.order;
                    return fail(message.str());
                }
            }
        }

        gpu.commitTrial(gpuTrial);
        const auto fullPostprocessorMs =
            gpuTrial.diagnostics.serializationPostprocessorMs +
            gpuTrial.diagnostics.energyCoenergyMs;
        auto residualOnly = gpu.evaluateTrial(
            femm::PhysicalResultLevel::ResidualOnly);
        if (!residualOnly.real ||
            residualOnly.real->physicalResultLevel !=
                femm::PhysicalResultLevel::ResidualOnly ||
            !residualOnly.real->airGaps.empty() ||
            residualOnly.real->magneticFieldEnergyJ != 0 ||
            residualOnly.real->magneticFieldCoenergyJ != 0 ||
            residualOnly.diagnostics.energyCoenergyMs != 0)
            return fail("residual-only evaluation performed accepted/full diagnostics");
        for (std::size_t i = 0; i < 3; ++i)
            if (!close(residualOnly.real->circuits[i].fluxLinkage,
                       gpuTrial.real->circuits[i].fluxLinkage, 1e-10, 1e-13))
                return fail("direct residual linkage changed validated semantics");
        const double residualPostprocessorMs =
            residualOnly.diagnostics.serializationPostprocessorMs +
            residualOnly.diagnostics.energyCoenergyMs;
        if (fullPostprocessorMs > 0 &&
            residualPostprocessorMs >= fullPostprocessorMs)
            return fail("residual-only evaluation did not skip full postprocessing");
        gpu.completeTrial(residualOnly,
                          femm::PhysicalResultLevel::AcceptedState);
        if (residualOnly.real->physicalResultLevel !=
                femm::PhysicalResultLevel::AcceptedState ||
            residualOnly.real->airGaps.empty() ||
            !close(residualOnly.real->torque, gpuTrial.real->torque,
                   1e-10, 1e-12))
            return fail("accepted-state promotion did not reproduce AGE torque");
        for (const auto &gapResult : residualOnly.real->airGaps)
            if (!gapResult.harmonics.empty())
                return fail("accepted-state promotion packaged diagnostic harmonics");
        gpu.commitTrial(residualOnly);
        const auto initialAngles = gpu.analysis().solveParameters().airGapPositions;
        const auto gap = initialAngles.begin()->first;
        const double originalInner = initialAngles.begin()->second.innerAngle;
        const double originalOuter = initialAngles.begin()->second.outerAngle;
        gpu.analysis().setAirGapAngle(gap, originalInner + 0.1, originalOuter);
        const auto rejected = gpu.evaluateTrial();
        if (!rejected.real || rejected.diagnostics.exactTopologyFallback ||
            rejected.diagnostics.symbolicAnalysisMs != 0 ||
            rejected.diagnostics.bucketConstructionMs != 0 ||
            rejected.diagnostics.bucketIdentity != gpuTrial.diagnostics.bucketIdentity)
            return fail("same-bucket transient trial did not use the immutable union");
        const double rejectedDifference = relativeVectorError(
            gpuTrial.real->nodal.magneticVectorPotential,
            rejected.real->nodal.magneticVectorPotential);
        if (rejectedDifference <= 1e-12)
            return fail("rejected angle trial did not produce a distinct trial field");
        gpu.rollbackToCommitted();
        gpu.analysis().setAirGapAngle(gap, originalInner + 360.0, originalOuter);
        if (!close(gpu.analysis().solveParameters().airGapPositions.at(gap).innerAngle,
                   originalInner, 0, 1e-12))
            return fail("+360 degree angle was not canonicalized");
        const auto restored = gpu.evaluateTrial();
        if (!restored.real)
            return fail("restored real solution payload missing");
        const double repeatError = relativeVectorError(
            gpuTrial.real->nodal.magneticVectorPotential,
            restored.real->nodal.magneticVectorPotential);
        if (repeatError > 1e-9)
            return fail("commit/rollback or modulo-angle repeatability regression: " +
                        std::to_string(repeatError));
        if (restored.diagnostics.bucketIdentity != gpuTrial.diagnostics.bucketIdentity ||
            restored.diagnostics.symbolicAnalysisMs != 0 ||
            restored.diagnostics.exactTopologyFallback) {
            std::ostringstream message;
            message << std::setprecision(17)
                    << "repeat evaluation did not reuse the canonical AGE bucket"
                    << " initial_bucket=" << gpuTrial.diagnostics.bucketIdentity
                    << " restored_bucket=" << restored.diagnostics.bucketIdentity
                    << " symbolic_ms=" << restored.diagnostics.symbolicAnalysisMs
                    << " exact_fallback="
                    << restored.diagnostics.exactTopologyFallback
                    << " bucket_reused=" << restored.diagnostics.bucketReused
                    << " symbolic_reused=" << restored.diagnostics.symbolicReused;
            return fail(message.str());
        }
        if (!close(restored.real->torque, gpuTrial.real->torque, 1e-7, 1e-12))
            return fail("rollback did not restore committed physical results");
        for (std::size_t i = 0; i < 3; ++i)
            if (!close(restored.real->circuits[i].fluxLinkage,
                       gpuTrial.real->circuits[i].fluxLinkage, 1e-7, 1e-12))
                return fail("rollback did not restore committed circuit results");
        if (relativeVectorError(rejected.real->nodal.magneticVectorPotential,
                                restored.real->nodal.magneticVectorPotential) <= 1e-12)
            return fail("rollback retained the rejected trial field");

        gpu.analysis().setAirGapAngle(gap, originalInner - 360.0, originalOuter);
        if (!close(gpu.analysis().solveParameters().airGapPositions.at(gap).innerAngle,
                   originalInner, 0, 1e-12))
            return fail("-360 degree angle was not canonicalized");
        const auto repeated = gpu.evaluateTrial();
        if (!repeated.real ||
            relativeVectorError(gpuTrial.real->nodal.magneticVectorPotential,
                                repeated.real->nodal.magneticVectorPotential) > 1e-9 ||
            repeated.diagnostics.symbolicAnalysisMs != 0 ||
            repeated.diagnostics.exactTopologyFallback ||
            !repeated.diagnostics.bucketReused)
            return fail("repeated canonical evaluation did not reuse committed resources");

        // Exercise the ESC-relevant case: one branch crosses exact zero while
        // the other two stay substantially loaded.  Check every returned
        // linkage, not only the linkage of the crossing branch.
        constexpr std::array<double, 5> crossingCurrents{
            -0.1, -0.01, 0.0, 0.01, 0.1};
        for (std::size_t crossing = 0; crossing < 3; ++crossing) {
            for (double backgroundSign : {1.0, -1.0}) {
                std::array<std::array<double, 3>, 5> linkages{};
                for (std::size_t sample = 0; sample < crossingCurrents.size();
                     ++sample) {
                    std::array<double, 3> currents{};
                    currents[crossing] = crossingCurrents[sample];
                    currents[(crossing + 1) % 3] = backgroundSign * 4.0;
                    currents[(crossing + 2) % 3] = backgroundSign * -2.0;
                    gpu.analysis().updateSolveParameters(
                        [&](femm::SolveParameters &parameters) {
                            for (std::size_t i = 0; i < currents.size(); ++i)
                                parameters.circuitConstraints[femm::CircuitId{i}] = {
                                    femm::CircuitConstraintKind::PrescribedCurrent,
                                    CComplex(currents[i], 0)};
                        });
                    const auto trial = gpu.evaluateTrial();
                    if (!trial.real || !trial.diagnostics.nonlinearWarmStartUsed)
                        return fail("loaded zero-crossing trial payload/warm start missing");
                    for (std::size_t i = 0; i < 3; ++i)
                        linkages[sample][i] = trial.real->circuits[i].fluxLinkage;
                }
                for (std::size_t linkage = 0; linkage < 3; ++linkage) {
                    const double negative = linkages[1][linkage];
                    const double zero = linkages[2][linkage];
                    const double positive = linkages[3][linkage];
                    const double midpointDefect =
                        std::abs(zero - 0.5 * (negative + positive));
                    const double localVariation = std::max(
                        {std::abs(positive - negative),
                         std::abs(zero - negative),
                         std::abs(positive - zero)});
                    if (midpointDefect > 5e-10 + 0.05 * localVariation)
                        return fail("loaded linkage is discontinuous through zero current");
                    const double leftSlope = (zero - negative) / 0.01;
                    const double rightSlope = (positive - zero) / 0.01;
                    const double slopeScale =
                        std::max(std::abs(leftSlope), std::abs(rightSlope));
                    if (std::abs(rightSlope - leftSlope) >
                        5e-8 + 0.15 * slopeScale)
                        return fail("loaded linkage slope is not smooth through zero current");
                }
            }
        }
        gpu.rollbackToCommitted();

        // Bound expensive device contexts while retaining host-side bucket
        // definitions. Start from a known bucket, reject a trial in its
        // neighbor, roll back, then commit monotonically through enough
        // buckets to force eviction. Nodal state lives at session scope, so
        // neither eviction nor bucket switching may discard the warm field.
        const std::array<double, 3> cacheCurrents{1.5, -0.75, 0.25};
        const auto setCacheState = [&](double relativeAngle) {
            gpu.analysis().updateSolveParameters(
                [&](femm::SolveParameters &parameters) {
                    parameters.airGapPositions[gap] = {
                        originalInner + relativeAngle, originalOuter};
                    for (std::size_t i = 0; i < cacheCurrents.size(); ++i)
                        parameters.circuitConstraints[femm::CircuitId{i}] = {
                            femm::CircuitConstraintKind::PrescribedCurrent,
                            CComplex(cacheCurrents[i], 0)};
                });
        };
        setCacheState(12.2); // bucket 10 for a 1.2-degree cache
        const auto cacheCommitted = gpu.evaluateTrial();
        if (!cacheCommitted.real ||
            !cacheCommitted.diagnostics.nonlinearWarmStartUsed)
            return fail("cache test committed-state warm start missing");
        gpu.commitTrial(cacheCommitted);

        setCacheState(13.4); // rejected trial across the 10 -> 11 boundary
        const auto crossBoundaryRejected = gpu.evaluateTrial();
        if (!crossBoundaryRejected.real ||
            !crossBoundaryRejected.diagnostics.nonlinearWarmStartUsed ||
            crossBoundaryRejected.diagnostics.bucketIdentity ==
                cacheCommitted.diagnostics.bucketIdentity ||
            gpu.residentBucketCount() > 2)
            return fail("cross-boundary rejected trial did not retain a bounded warm cache");
        gpu.rollbackToCommitted();
        if (gpu.residentBucketCount() > 2)
            return fail("rollback exceeded the cuDSS bucket cache capacity");
        setCacheState(12.2);
        const auto cacheRestored = gpu.evaluateTrial();
        if (!cacheRestored.real ||
            cacheRestored.diagnostics.bucketIdentity !=
                cacheCommitted.diagnostics.bucketIdentity ||
            !cacheRestored.diagnostics.bucketReused ||
            !cacheRestored.diagnostics.symbolicReused ||
            !cacheRestored.diagnostics.nonlinearWarmStartUsed ||
            relativeVectorError(
                cacheCommitted.real->nodal.magneticVectorPotential,
                cacheRestored.real->nodal.magneticVectorPotential) > 1e-9)
            return fail("cross-boundary rollback did not restore the committed warm state");
        gpu.commitTrial(cacheRestored);

        const std::size_t evictionsBefore = gpu.bucketEvictionCount();
        femm::TrialSolution lastCommitted = cacheRestored;
        for (double relativeAngle : {13.4, 14.6, 15.8, 17.0}) {
            setCacheState(relativeAngle);
            auto next = gpu.evaluateTrial();
            if (!next.real || !next.diagnostics.nonlinearWarmStartUsed ||
                next.diagnostics.residentBucketCount > 2 ||
                gpu.residentBucketCount() > 2)
                return fail("monotonic bucket crossing lost warm start or cache bound");
            gpu.commitTrial(next);
            lastCommitted = std::move(next);
        }
        if (gpu.bucketEvictionCount() <= evictionsBefore ||
            gpu.residentBucketCount() != 2 ||
            gpu.bucketDefinitionCount() < 5)
            return fail("bounded cache did not evict device contexts while retaining definitions");

        // Revisit an evicted bucket. It should reanalyze its retained host
        // definition, while the last committed bucket stays resident for a
        // subsequent rollback.
        setCacheState(12.2);
        const auto recreated = gpu.evaluateTrial();
        if (!recreated.real || recreated.diagnostics.bucketReused ||
            recreated.diagnostics.symbolicAnalysisMs <= 0 ||
            !recreated.diagnostics.nonlinearWarmStartUsed ||
            gpu.residentBucketCount() != 2)
            return fail("evicted bucket was not safely recreated from host metadata");
        gpu.rollbackToCommitted();
        setCacheState(17.0);
        const auto finalRestored = gpu.evaluateTrial();
        if (!finalRestored.real || !finalRestored.diagnostics.bucketReused ||
            !finalRestored.diagnostics.symbolicReused ||
            relativeVectorError(
                lastCommitted.real->nodal.magneticVectorPotential,
                finalRestored.real->nodal.magneticVectorPotential) > 1e-9)
            return fail("rollback failed after recreating an evicted trial bucket");

        if (gpu.analysis().meshGenerationCount() != 1)
            return fail("persistent evaluations unexpectedly regenerated the mesh");
        if (gpu.topologyImportCount() != 1 || gpu.orderingCount() != 1 ||
            gpu.meshFileReadCount() != 0 || gpu.meshFileWriteCount() != 0)
            return fail("persistent mesh/order state or in-memory I/O contract regressed");
    } catch (const std::exception &error) {
        return fail(error.what());
    }
    return 0;
}
