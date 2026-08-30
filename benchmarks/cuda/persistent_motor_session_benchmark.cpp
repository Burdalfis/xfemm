#include "PersistentMotorSession.h"

#include "FemmReader.h"
#include "linsolve/CudssLinearSystemBackend.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::shared_ptr<femm::FemmProblem> loadProblem(const std::string &path)
{
    auto problem = std::make_shared<femm::FemmProblem>(
        femm::FileType::MagneticsFile);
    std::ostringstream errors;
    femm::MagneticsReader reader(problem, errors);
    if (reader.parse(path) != femm::F_FILE_OK)
        throw std::runtime_error("could not parse motor fixture: " + errors.str());
    return problem;
}

double sessionOtherMs(const femm::EvaluationDiagnostics &d)
{
    return std::max(0.0, d.sessionSynchronizationMs - d.airGapUpdateMs -
                           d.bucketDefinitionConstructionMs);
}

double accountedMs(const femm::EvaluationDiagnostics &d)
{
    return sessionOtherMs(d) + d.modelPreparationMs + d.bucketLookupMs +
           d.bucketSwitchMs + d.bucketConstructionMs + d.symbolicAnalysisMs +
           d.airGapUpdateMs + d.nonlinearMaterialEvaluationMs +
           d.numericMatrixAssemblyMs + d.sparsePackingMs + d.hostToDeviceMs +
           d.numericFactorizationMs + d.linearSolveMs + d.deviceToHostMs +
           d.residualEvaluationMs + d.nonlinearBookkeepingMs +
           d.resultPackagingMs + d.serializationPostprocessorMs +
           d.fluxLinkageMs + d.torqueMs +
           d.airGapHarmonicPackagingMs + d.energyCoenergyMs;
}

void printHeader()
{
    std::cout
        << "PROFILE_HEADER phase,index,angle_deg,total_ms,session_other_ms,"
           "model_prepare_ms,bucket_lookup_ms,bucket_switch_ms,"
           "bucket_definition_ms,bucket_resource_ms,symbolic_ms,age_update_ms,"
           "material_ms,matrix_assembly_ms,sparse_pack_scatter_ms,h2d_ms,"
           "age_matrix_ms,element_matrix_and_rhs_ms,explicit_rhs_ms,"
           "boundary_constraints_ms,"
           "device_clear_ms,device_material_ms,device_element_ms,"
           "device_scatter_ms,device_age_upload_ms,device_constraint_ms,"
           "factor_ms,solve_ms,d2h_ms,residual_ms,nonlinear_bookkeeping_ms,"
           "result_packaging_ms,postprocessor_update_ms,flux_ms,torque_ms,"
           "age_harmonic_packaging_ms,energy_coenergy_ms,"
           "unaccounted_ms,newton_iterations,matrix_nonzeros,factor_nonzeros,"
           "permanent_device_bytes,peak_device_bytes,h2d_bytes,relative_residual,"
           "bucket_reused,symbolic_reused,exact_fallback,device_assembly,"
           "parity_max_abs_entry,parity_max_rel_entry,parity_max_abs_rhs,"
           "parity_symmetry,bucket_identity\n";
}

void printProfile(const std::string &phase, std::size_t index, double angle,
                  const femm::EvaluationDiagnostics &d)
{
    const double bucketResource =
        std::max(0.0, d.bucketConstructionMs - d.bucketDefinitionConstructionMs);
    const double unaccounted = std::max(0.0, d.totalEvaluateMs - accountedMs(d));
    std::cout << std::setprecision(9) << "PROFILE " << phase << ',' << index << ','
              << angle << ',' << d.totalEvaluateMs << ',' << sessionOtherMs(d) << ','
              << d.modelPreparationMs << ',' << d.bucketLookupMs << ','
              << d.bucketSwitchMs << ',' << d.bucketDefinitionConstructionMs << ','
              << bucketResource << ',' << d.symbolicAnalysisMs << ','
              << d.airGapUpdateMs << ',' << d.nonlinearMaterialEvaluationMs << ','
              << d.numericMatrixAssemblyMs << ',' << d.sparsePackingMs << ','
              << d.hostToDeviceMs << ',' << d.airGapMatrixAssemblyMs << ','
              << d.elementMatrixAssemblyMs << ',' << d.rhsConstructionMs << ','
              << d.boundaryConditionApplicationMs << ','
              << d.deviceAssemblyClearMs << ',' << d.deviceMaterialMs << ','
              << d.deviceElementMs << ',' << d.deviceScatterMs << ','
              << d.deviceAgeUploadMs << ',' << d.deviceConstraintMs << ','
              << d.numericFactorizationMs << ','
              << d.linearSolveMs << ',' << d.deviceToHostMs << ','
              << d.residualEvaluationMs << ',' << d.nonlinearBookkeepingMs << ','
              << d.resultPackagingMs << ',' << d.serializationPostprocessorMs << ','
              << d.fluxLinkageMs << ',' << d.torqueMs << ','
              << d.airGapHarmonicPackagingMs << ',' << d.energyCoenergyMs << ','
              << unaccounted << ','
              << d.nonlinearIterations << ',' << d.matrixNonzeros << ','
              << d.factorNonzeros << ','
              << d.permanentDeviceBytes << ',' << d.peakDeviceBytes << ','
              << d.hostToDeviceBytes << ',' << d.relativeResidual << ','
              << d.bucketReused << ',' << d.symbolicReused << ','
              << d.exactTopologyFallback << ',' << d.deviceAssemblyUsed << ','
              << d.assemblyParityMaxAbsoluteEntry << ','
              << d.assemblyParityMaxRelativeEntry << ','
              << d.assemblyParityMaxAbsoluteRhs << ','
              << d.assemblyParitySymmetryDifference << ','
              << d.bucketIdentity << '\n';
}

femm::TrialSolution evaluateOuterAngle(
    femm::PersistentMotorSession &session, double outerAngle,
    const std::array<double, 3> &currents,
    femm::PhysicalResultLevel level = femm::PhysicalResultLevel::FullDiagnostics)
{
    session.analysis().updateSolveParameters([&](femm::SolveParameters &parameters) {
        for (std::size_t i = 0; i < currents.size(); ++i)
            parameters.circuitConstraints[femm::CircuitId{i}] = {
                femm::CircuitConstraintKind::PrescribedCurrent,
                CComplex(currents[i], 0)};
        for (auto &gap : parameters.airGapPositions)
            gap.second.outerAngle = outerAngle;
    });
    return session.evaluateTrial(level);
}

femm::PhysicalResultLevel parseResultLevel(const std::string &name)
{
    if (name == "residual") return femm::PhysicalResultLevel::ResidualOnly;
    if (name == "accepted") return femm::PhysicalResultLevel::AcceptedState;
    if (name == "full") return femm::PhysicalResultLevel::FullDiagnostics;
    throw std::invalid_argument("evaluation level must be residual, accepted, or full");
}

femm::PlanarAssemblyBackend parseAssemblyBackend(const std::string &name)
{
    if (name == "host") return femm::PlanarAssemblyBackend::Host;
    if (name == "cuda") return femm::PlanarAssemblyBackend::CudaAtomic;
    if (name == "cuda-deterministic")
        return femm::PlanarAssemblyBackend::CudaDeterministic;
    throw std::invalid_argument(
        "assembly backend must be host, cuda, or cuda-deterministic");
}

void printPhysical(const std::string &phase, const femm::TrialSolution &result)
{
    std::cout << std::setprecision(17) << "PHYSICAL phase=" << phase
              << " nodes=" << result.real->nodal.magneticVectorPotential.size()
              << " torque=" << result.real->torque;
    for (const auto &circuit : result.real->circuits)
        std::cout << " flux=" << circuit.fluxLinkage;
    std::cout << " air_gaps=" << result.real->airGaps.size();
    for (const auto &gap : result.real->airGaps)
        std::cout << " age_torque=" << gap.torque
                  << " harmonics=" << gap.harmonics.size();
    std::cout << '\n';
}

template<typename Member>
double mean(const std::vector<femm::EvaluationDiagnostics> &samples, Member member)
{
    double total = 0;
    for (const auto &sample : samples) total += sample.*member;
    return samples.empty() ? 0 : total / samples.size();
}

void printHotSummary(const std::vector<femm::EvaluationDiagnostics> &samples)
{
    if (samples.empty()) return;
    double totalMin = samples.front().totalEvaluateMs;
    double totalMax = totalMin;
    double unaccounted = 0;
    double sessionOther = 0;
    double bucketResource = 0;
    for (const auto &sample : samples) {
        totalMin = std::min(totalMin, sample.totalEvaluateMs);
        totalMax = std::max(totalMax, sample.totalEvaluateMs);
        unaccounted += std::max(0.0, sample.totalEvaluateMs - accountedMs(sample));
        sessionOther += sessionOtherMs(sample);
        bucketResource += std::max(
            0.0, sample.bucketConstructionMs - sample.bucketDefinitionConstructionMs);
    }
    const double count = static_cast<double>(samples.size());
    std::cout << std::setprecision(9)
              << "PROFILE_SUMMARY phase=transient_hot samples=" << samples.size()
              << " total_mean_ms=" << mean(samples, &femm::EvaluationDiagnostics::totalEvaluateMs)
              << " total_min_ms=" << totalMin
              << " total_max_ms=" << totalMax
              << " session_other_mean_ms=" << sessionOther / count
              << " model_prepare_mean_ms="
              << mean(samples, &femm::EvaluationDiagnostics::modelPreparationMs)
              << " bucket_lookup_mean_ms="
              << mean(samples, &femm::EvaluationDiagnostics::bucketLookupMs)
              << " bucket_switch_mean_ms="
              << mean(samples, &femm::EvaluationDiagnostics::bucketSwitchMs)
              << " bucket_definition_mean_ms="
              << mean(samples, &femm::EvaluationDiagnostics::bucketDefinitionConstructionMs)
              << " bucket_resource_mean_ms=" << bucketResource / count
              << " symbolic_mean_ms="
              << mean(samples, &femm::EvaluationDiagnostics::symbolicAnalysisMs)
              << " age_update_mean_ms="
              << mean(samples, &femm::EvaluationDiagnostics::airGapUpdateMs)
              << " material_mean_ms="
              << mean(samples, &femm::EvaluationDiagnostics::nonlinearMaterialEvaluationMs)
              << " matrix_assembly_mean_ms="
              << mean(samples, &femm::EvaluationDiagnostics::numericMatrixAssemblyMs)
              << " age_matrix_mean_ms="
              << mean(samples, &femm::EvaluationDiagnostics::airGapMatrixAssemblyMs)
              << " element_matrix_and_rhs_mean_ms="
              << mean(samples, &femm::EvaluationDiagnostics::elementMatrixAssemblyMs)
              << " explicit_rhs_mean_ms="
              << mean(samples, &femm::EvaluationDiagnostics::rhsConstructionMs)
              << " boundary_constraints_mean_ms="
              << mean(samples, &femm::EvaluationDiagnostics::boundaryConditionApplicationMs)
              << " device_clear_mean_ms="
              << mean(samples, &femm::EvaluationDiagnostics::deviceAssemblyClearMs)
              << " device_material_mean_ms="
              << mean(samples, &femm::EvaluationDiagnostics::deviceMaterialMs)
              << " device_element_mean_ms="
              << mean(samples, &femm::EvaluationDiagnostics::deviceElementMs)
              << " device_scatter_mean_ms="
              << mean(samples, &femm::EvaluationDiagnostics::deviceScatterMs)
              << " device_age_upload_mean_ms="
              << mean(samples, &femm::EvaluationDiagnostics::deviceAgeUploadMs)
              << " device_constraint_mean_ms="
              << mean(samples, &femm::EvaluationDiagnostics::deviceConstraintMs)
              << " sparse_pack_scatter_mean_ms="
              << mean(samples, &femm::EvaluationDiagnostics::sparsePackingMs)
              << " h2d_mean_ms="
              << mean(samples, &femm::EvaluationDiagnostics::hostToDeviceMs)
              << " factor_mean_ms="
              << mean(samples, &femm::EvaluationDiagnostics::numericFactorizationMs)
              << " solve_mean_ms="
              << mean(samples, &femm::EvaluationDiagnostics::linearSolveMs)
              << " d2h_mean_ms="
              << mean(samples, &femm::EvaluationDiagnostics::deviceToHostMs)
              << " residual_mean_ms="
              << mean(samples, &femm::EvaluationDiagnostics::residualEvaluationMs)
              << " nonlinear_bookkeeping_mean_ms="
              << mean(samples, &femm::EvaluationDiagnostics::nonlinearBookkeepingMs)
              << " result_packaging_mean_ms="
              << mean(samples, &femm::EvaluationDiagnostics::resultPackagingMs)
              << " postprocessor_update_mean_ms="
              << mean(samples, &femm::EvaluationDiagnostics::serializationPostprocessorMs)
              << " flux_mean_ms="
              << mean(samples, &femm::EvaluationDiagnostics::fluxLinkageMs)
              << " torque_mean_ms="
              << mean(samples, &femm::EvaluationDiagnostics::torqueMs)
              << " age_harmonic_packaging_mean_ms="
              << mean(samples, &femm::EvaluationDiagnostics::airGapHarmonicPackagingMs)
              << " energy_coenergy_mean_ms="
              << mean(samples, &femm::EvaluationDiagnostics::energyCoenergyMs)
              << " unaccounted_mean_ms=" << unaccounted / count
              << " newton_iterations_mean="
              << mean(samples, &femm::EvaluationDiagnostics::nonlinearIterations)
              << '\n';
}

} // namespace

int main(int argc, char **argv)
{
    if (argc < 2 || argc > 5) {
        std::cerr << "expected a sliding-AGE motor .fem path, optional bucket "
                     "width, optional residual|accepted|full level, and optional "
                     "host|cuda|cuda-deterministic assembly backend\n";
        return 2;
    }
    try {
        const bool legacyReference = argc == 3 && std::string(argv[2]) == "legacy";
        const double bucketWidth = argc >= 3 && !legacyReference
            ? std::stod(argv[2]) : 1.2;
        const auto resultLevel = argc >= 4
            ? parseResultLevel(argv[3])
            : femm::PhysicalResultLevel::FullDiagnostics;
        const auto assemblyBackend = argc == 5
            ? parseAssemblyBackend(argv[4])
            : femm::PlanarAssemblyBackend::Host;
        if (legacyReference) {
            femm::PersistentMotorSession reference{
                femm::ModelDefinition(loadProblem(argv[1]))};
            const auto &parameters = reference.analysis().solveParameters();
            const double angle =
                parameters.airGapPositions.begin()->second.outerAngle;
            std::array<double, 3> currents{};
            for (std::size_t i = 0; i < currents.size(); ++i)
                currents[i] = parameters.circuitConstraints.at({i}).value.re;
            auto result = evaluateOuterAngle(reference, angle, currents);
            std::cout << std::setprecision(17)
                      << "LEGACY_REFERENCE nodes="
                      << result.real->nodal.magneticVectorPotential.size()
                      << " newton_iterations=" << result.diagnostics.nonlinearIterations
                      << " total_ms=" << result.diagnostics.totalEvaluateMs
                      << " torque=" << result.real->torque;
            for (const auto &circuit : result.real->circuits)
                std::cout << " flux=" << circuit.fluxLinkage;
            std::cout << '\n';
            return 0;
        }
        if (!femm::cudssDeviceAvailable()) {
            std::cerr << "No CUDA device\n";
            return 77;
        }
        femm::PersistentMotorSession session(
            femm::ModelDefinition(loadProblem(argv[1])),
            femm::CudssSessionOptions{bucketWidth, false, false, 2,
                                      assemblyBackend});
        if (std::abs(session.analysis().model().problem().Precision - 1e-8) > 1e-15)
            throw std::runtime_error("full-size production profile requires 1e-8 precision");
        if (session.analysis().solveParameters().airGapPositions.empty() ||
            session.analysis().solveParameters().circuitConstraints.size() != 3)
            throw std::runtime_error("benchmark fixture must have one AGE and three circuits");

        const auto &parameters = session.analysis().solveParameters();
        const double originalAngle =
            parameters.airGapPositions.begin()->second.outerAngle;
        std::array<double, 3> baseCurrents{};
        for (std::size_t i = 0; i < baseCurrents.size(); ++i)
            baseCurrents[i] = parameters.circuitConstraints.at({i}).value.re;

        session.initialize();
        printHeader();
        printProfile("initialization", 0, originalAngle,
                     session.initializationDiagnostics());

        auto baseline = evaluateOuterAngle(
            session, originalAngle, baseCurrents, resultLevel);
        const std::size_t nodes =
            baseline.real->nodal.magneticVectorPotential.size();
        const std::size_t storedNonzeros = baseline.diagnostics.matrixNonzeros;
        const std::size_t fullSymmetricNonzeros = 2 * storedNonzeros - nodes;
        std::cout << "WORKLOAD nodes=" << nodes
                  << " stored_triangle_nnz=" << storedNonzeros
                  << " full_symmetric_nnz=" << fullSymmetricNonzeros << '\n';
        if (nodes < 145000 || nodes > 155000 ||
            fullSymmetricNonzeros < 1050000 || fullSymmetricNonzeros > 1100000)
            throw std::runtime_error(
                "persistent benchmark did not produce the expected full-size motor system");
        printProfile("primed_hot", 0, originalAngle, baseline.diagnostics);
        printPhysical("primed_hot", baseline);
        session.commitTrial(baseline);

        std::vector<femm::EvaluationDiagnostics> hotSamples;
        std::array<double, 3> currents = baseCurrents;
        constexpr double twoPi = 6.28318530717958647692;
        for (std::size_t i = 0; i < 12; ++i) {
            const double angle = originalAngle + 0.01 * static_cast<double>(i + 1);
            const double phase = 0.08 * static_cast<double>(i + 1);
            for (std::size_t c = 0; c < currents.size(); ++c)
                currents[c] = baseCurrents[c] +
                    0.25 * std::sin(phase - twoPi * static_cast<double>(c) / 3.0);
            auto trial = evaluateOuterAngle(
                session, angle, currents, resultLevel);
            printProfile("transient_hot", i, angle, trial.diagnostics);
            hotSamples.push_back(trial.diagnostics);
            session.commitTrial(trial);
        }
        printHotSummary(hotSamples);

        const double switchedAngle = originalAngle + 0.25;
        auto switched = evaluateOuterAngle(
            session, switchedAngle, currents, resultLevel);
        printProfile("new_bucket", 0, switchedAngle, switched.diagnostics);
        session.commitTrial(switched);

        const double returnAngle = originalAngle + 0.08;
        auto returned = evaluateOuterAngle(
            session, returnAngle, currents, resultLevel);
        printProfile("cached_bucket_return", 0, returnAngle, returned.diagnostics);
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
