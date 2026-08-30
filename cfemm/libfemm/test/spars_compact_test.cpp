#include "spars.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>

namespace {

bool close(double actual, double expected, double tolerance = 1.e-12)
{
    return std::abs(actual - expected) <= tolerance;
}

bool checkVector(const char *description, const std::array<double, 3> &actual,
                 const std::array<double, 3> &expected)
{
    for (std::size_t i = 0; i < actual.size(); ++i) {
        if (!close(actual[i], expected[i])) {
            std::cerr << description << '[' << i << "]: expected " << expected[i]
                      << ", got " << actual[i] << '\n';
            return false;
        }
    }
    return true;
}

} // namespace

int main()
{
    CBigLinProb problem;
    if (!problem.Create(3, 0))
        return 1;

    // Symmetric matrix [[4,1,0],[1,3,2],[0,2,5]].
    problem.Put(4, 0, 0);
    problem.Put(1, 0, 1);
    problem.Put(3, 1, 1);
    problem.Put(2, 1, 2);
    problem.Put(5, 2, 2);

    const std::array<double, 3> input{{1, 2, 3}};
    std::array<double, 3> product{};
    problem.MultA(input.data(), product.data());
    bool ok = checkVector("initial product", product, {{6, 13, 19}});

    // Verify that mutation invalidates and refreshes the packed values.
    problem.Wipe();
    problem.AddTo(2, 0, 0);
    problem.AddTo(1, 1, 1);
    problem.AddTo(4, 2, 2);
    problem.MultA(input.data(), product.data());
    ok &= checkVector("rebuilt product", product, {{2, 2, 12}});

    // The persistent element path must be numerically identical to six
    // ordinary sparse additions and must retain valid destinations after a
    // numerical wipe.
    CBigLinProb direct;
    CBigLinProb reference;
    if (!direct.Create(3, 0) || !reference.Create(3, 0))
        return 1;
    const int elementNodes[3] = {0, 1, 2};
    const double firstElement[6] = {4., 1., 0., 3., 2., 5.};
    const double firstRhs[3] = {1.25, -2.5, .75};
    direct.AddSymmetric3x3(0, elementNodes, firstElement);
    std::size_t entry = 0;
    for (int row = 0; row < 3; ++row)
        for (int column = row; column < 3; ++column)
            reference.AddTo(firstElement[entry++], row, column);
    for (int row = 0; row < 3; ++row) {
        direct.b[row] += firstRhs[row];
        reference.b[row] += firstRhs[row];
    }
    std::vector<std::int32_t> directRows, directColumns;
    std::vector<std::int32_t> referenceRows, referenceColumns;
    std::vector<double> directValues, referenceValues;
    direct.copyUpperCsr(directRows, directColumns, directValues);
    reference.copyUpperCsr(referenceRows, referenceColumns, referenceValues);
    if (directRows != referenceRows || directColumns != referenceColumns ||
        directValues != referenceValues) {
        std::cerr << "direct element assembly differs from scalar reference\n";
        ok = false;
    }
    direct.Wipe();
    reference.Wipe();
    const double secondElement[6] = {2., -1., .5, 7., 3., 4.};
    const double secondRhs[3] = {-3., .125, 8.};
    direct.AddSymmetric3x3(0, elementNodes, secondElement);
    entry = 0;
    for (int row = 0; row < 3; ++row)
        for (int column = row; column < 3; ++column)
            reference.AddTo(secondElement[entry++], row, column);
    for (int row = 0; row < 3; ++row) {
        direct.b[row] += secondRhs[row];
        reference.b[row] += secondRhs[row];
    }
    direct.copyUpperCsr(directRows, directColumns, directValues);
    reference.copyUpperCsr(referenceRows, referenceColumns, referenceValues);
    if (directRows != referenceRows || directColumns != referenceColumns ||
        directValues != referenceValues) {
        std::cerr << "retained element destinations differ after wipe\n";
        ok = false;
    }
    double maximumAbsoluteEntryDifference = 0.;
    double maximumRelativeEntryDifference = 0.;
    for (std::size_t i = 0; i < directValues.size(); ++i) {
        const double difference = std::abs(directValues[i] - referenceValues[i]);
        maximumAbsoluteEntryDifference = std::max(
            maximumAbsoluteEntryDifference, difference);
        maximumRelativeEntryDifference = std::max(
            maximumRelativeEntryDifference,
            difference / std::max(1.e-300, std::abs(referenceValues[i])));
    }
    double maximumRhsDifference = 0.;
    double maximumSymmetryDifference = 0.;
    for (int row = 0; row < 3; ++row) {
        maximumRhsDifference = std::max(
            maximumRhsDifference, std::abs(direct.b[row] - reference.b[row]));
        for (int column = 0; column < 3; ++column)
            maximumSymmetryDifference = std::max(
                maximumSymmetryDifference,
                std::abs(direct.Get(row, column) - direct.Get(column, row)));
    }
    std::cout << "spars_direct_assembly_parity max_abs_entry_difference="
              << maximumAbsoluteEntryDifference
              << " max_relative_entry_difference="
              << maximumRelativeEntryDifference
              << " max_rhs_difference=" << maximumRhsDifference
              << " max_symmetry_difference=" << maximumSymmetryDifference
              << '\n';

    // Restore the SPD matrix and exercise the packed SSOR/PCG path.
    problem.Wipe();
    problem.AddTo(4, 0, 0);
    problem.AddTo(1, 0, 1);
    problem.AddTo(3, 1, 1);
    problem.AddTo(2, 1, 2);
    problem.AddTo(5, 2, 2);
    problem.b[0] = 6;
    problem.b[1] = 13;
    problem.b[2] = 19;
    problem.Precision = 1.e-12;
    ok &= problem.PCGSolve(0);
    if (problem.lastIterations() <= 0 ||
        problem.lastRelativeResidual() > problem.Precision) {
        std::cerr << "invalid PCG diagnostics: iterations="
                  << problem.lastIterations() << ", relative residual="
                  << problem.lastRelativeResidual() << '\n';
        ok = false;
    }
    for (std::size_t i = 0; i < input.size(); ++i) {
        if (!close(problem.V[i], input[i], 1.e-10)) {
            std::cerr << "solution[" << i << "]: expected " << input[i]
                      << ", got " << problem.V[i] << '\n';
            ok = false;
        }
    }

    // Exercise the escape representation at the reserved 16-bit sentinel.
    const char *columnIndex = std::getenv("XFEMM_PCG_COLUMN_INDEX");
    if (columnIndex != nullptr &&
        (std::strcmp(columnIndex, "mixed16") == 0 ||
         std::strcmp(columnIndex, "row16") == 0)) {
        constexpr int last = 65535;
        CBigLinProb wideProblem;
        if (!wideProblem.Create(last + 1, 0))
            return 1;
        wideProblem.Put(2., 0, last);
        std::vector<double> wideInput(static_cast<std::size_t>(last) + 1, 0.);
        std::vector<double> wideProduct(wideInput.size(), 0.);
        wideInput[0] = 3.;
        wideInput[static_cast<std::size_t>(last)] = 5.;
        wideProblem.MultA(wideInput.data(), wideProduct.data());
        if (!close(wideProduct[0], 10.) ||
            !close(wideProduct[static_cast<std::size_t>(last)], 6.)) {
            std::cerr << "16-bit wide-column product: expected [10, 6], got ["
                      << wideProduct[0] << ", "
                      << wideProduct[static_cast<std::size_t>(last)] << "]\n";
            ok = false;
        }
    }

    return ok ? 0 : 1;
}
