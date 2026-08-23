#include "spars.h"

#include <array>
#include <cmath>
#include <iostream>

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
    for (std::size_t i = 0; i < input.size(); ++i) {
        if (!close(problem.V[i], input[i], 1.e-10)) {
            std::cerr << "solution[" << i << "]: expected " << input[i]
                      << ", got " << problem.V[i] << '\n';
            ok = false;
        }
    }

    return ok ? 0 : 1;
}
