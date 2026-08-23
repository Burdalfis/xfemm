#include "FemmProblem.h"

#include <cmath>
#include <iostream>
#include <memory>

int main()
{
    femm::FemmProblem problem(femm::FileType::MagneticsFile);
    problem.nodelist.push_back(std::make_unique<femm::CNode>(3.05, 0.));
    problem.nodelist.push_back(std::make_unique<femm::CNode>(-3.05, 0.));

    femm::CArcSegment semicircle;
    semicircle.n0 = 0;
    semicircle.n1 = 1;
    semicircle.ArcLength = 180.;

    CComplex center;
    double radius = 0.;
    problem.getCircle(semicircle, center, radius);

    if (center.re != 0. || center.im != 0.) {
        std::cerr << "semicircle center: expected (0,0), got ("
                  << center.re << ',' << center.im << ")\n";
        return 1;
    }
    if (std::abs(radius - 3.05) > 1.e-15) {
        std::cerr << "semicircle radius: expected 3.05, got " << radius << '\n';
        return 1;
    }

    return 0;
}
