#include "CMaterialProp.h"
#include "linsolve/CudaPlanarAssembly.h"
#include "linsolve/CudssLinearSystemBackend.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

int main()
{
    if (!femm::cudssDeviceAvailable()) return 77;

    femm::CMSolverMaterialProp host;
    host.Bdata = {0.0, 0.35, 0.8, 1.2, 1.55, 1.9, 2.2};
    host.Hdata = {0.0, 28.0, 70.0, 145.0, 420.0, 3200.0, 18000.0};
    host.slope = {70.0, 82.0, 130.0, 420.0, 3200.0, 22000.0, 52000.0};
    host.BHpoints = static_cast<int>(host.Bdata.size());

    femm::PlanarAssemblyMaterial material;
    material.bhCount = host.BHpoints;
    std::vector<double> field;
    std::vector<double> slope;
    for (const auto &value : host.Hdata) field.push_back(value.re);
    for (const auto &value : host.slope) slope.push_back(value.re);
    const std::vector<double> samples = {
        0.0, 1e-12, -0.05, 0.35, -0.6, 0.8, 1.0, 1.2,
        -1.4, 1.55, 1.75, 1.9, -2.2, 2.5};
    const auto device = femm::evaluateCudaPlanarBh(
        material, host.Bdata, field, slope, samples);
    double maxV = 0.;
    double maxDv = 0.;
    for (std::size_t i = 0; i < samples.size(); ++i) {
        double expectedV = 0.;
        double expectedDv = 0.;
        host.GetBHProps(samples[i], expectedV, expectedDv);
        maxV = std::max(maxV, std::abs(expectedV - device[i].reluctivity));
        maxDv = std::max(maxDv,
            std::abs(expectedDv - device[i].differentialReluctivity));
    }
    std::cout << "cuda_planar_material_parity max_reluctivity_abs=" << maxV
              << " max_differential_abs=" << maxDv << '\n';
    if (maxV > 1e-12 || maxDv > 1e-12)
        throw std::runtime_error("CUDA and host B-H interpolation differ");
    return 0;
}
