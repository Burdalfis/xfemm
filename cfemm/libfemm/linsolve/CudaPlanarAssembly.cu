#include "CudaPlanarAssembly.h"

#include <cuda_runtime.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace femm {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kEquationScale = kPi * 4.e-5;
constexpr double kMu0 = kPi * 4.e-7;

void checkCuda(cudaError_t status, const char *operation)
{
    if (status != cudaSuccess)
        throw std::runtime_error(std::string(operation) + ": " +
                                 cudaGetErrorString(status));
}

#define XFEMM_CUDA_ASSEMBLY_CHECK(expr) checkCuda((expr), #expr)

template <typename T>
class Buffer
{
public:
    Buffer() = default;
    explicit Buffer(std::size_t count) { resize(count); }
    ~Buffer() { if (m_data) cudaFree(m_data); }
    Buffer(const Buffer &) = delete;
    Buffer &operator=(const Buffer &) = delete;

    void resize(std::size_t count)
    {
        if (count == m_count) return;
        if (m_data) XFEMM_CUDA_ASSEMBLY_CHECK(cudaFree(m_data));
        m_data = nullptr;
        m_count = count;
        if (count)
            XFEMM_CUDA_ASSEMBLY_CHECK(cudaMalloc(&m_data, count * sizeof(T)));
    }
    T *get() const { return m_data; }
    std::size_t size() const { return m_count; }
    std::size_t bytes() const { return m_count * sizeof(T); }

private:
    T *m_data = nullptr;
    std::size_t m_count = 0;
};

class Event
{
public:
    Event() { XFEMM_CUDA_ASSEMBLY_CHECK(cudaEventCreate(&m_event)); }
    ~Event() { if (m_event) cudaEventDestroy(m_event); }
    Event(const Event &) = delete;
    Event &operator=(const Event &) = delete;
    cudaEvent_t get() const { return m_event; }

private:
    cudaEvent_t m_event = nullptr;
};

double elapsedDevice(const Event &begin, const Event &end)
{
    float milliseconds = 0;
    XFEMM_CUDA_ASSEMBLY_CHECK(
        cudaEventElapsedTime(&milliseconds, begin.get(), end.get()));
    return milliseconds;
}

std::int32_t findLowerEntry(const std::vector<std::int32_t> &rows,
                            const std::vector<std::int32_t> &columns,
                            std::int32_t first, std::int32_t second)
{
    const std::int32_t row = std::max(first, second);
    const std::int32_t column = std::min(first, second);
    auto begin = columns.begin() + rows[static_cast<std::size_t>(row)];
    auto end = columns.begin() + rows[static_cast<std::size_t>(row) + 1];
    auto found = std::lower_bound(begin, end, column);
    if (found == end || *found != column)
        throw std::logic_error("CUDA planar assembly destination is absent from bucket CSR");
    return static_cast<std::int32_t>(found - columns.begin());
}

struct DeviceElement
{
    std::int32_t nodes[3];
    double p[3];
    double q[3];
    double mx[9];
    double my[9];
    double mxy[9];
    double fixedMatrix[9];
    double fixedRhs[3];
    double area;
    std::int32_t material;
    std::int32_t circuit;
};

struct DeviceMaterial
{
    double muX;
    double muY;
    double fill;
    double conductivity;
    std::int32_t laminationType;
    std::int32_t bhOffset;
    std::int32_t bhCount;
};

struct DeviceAdd
{
    std::int32_t destination;
    double value;
};

__device__ void bhProperties(const DeviceMaterial &material, double B,
                             const double *flux, const double *field,
                             const double *slope, double &v, double &dv)
{
    const double b = fabs(B);
    if (material.bhCount == 0) {
        v = material.muX;
        dv = 0.;
        return;
    }
    const int offset = material.bhOffset;
    if (b == 0.) {
        v = slope[offset];
        dv = 0.;
        return;
    }
    const int last = offset + material.bhCount - 1;
    if (b > flux[last]) {
        const double h = field[last] + slope[last] * (b - flux[last]);
        const double dh = slope[last];
        v = h / b;
        dv = 0.5 * (dh / (b * b) - h / (b * b * b));
        return;
    }
    for (int local = 0; local < material.bhCount - 1; ++local) {
        const int i = offset + local;
        if (b < flux[i] || b > flux[i + 1]) continue;
        const double length = flux[i + 1] - flux[i];
        const double z = (b - flux[i]) / length;
        const double z2 = z * z;
        const double h = (1. - 3. * z2 + 2. * z2 * z) * field[i] +
            z * (1. - 2. * z + z2) * length * slope[i] +
            z2 * (3. - 2. * z) * field[i + 1] +
            z2 * (z - 1.) * length * slope[i + 1];
        const double dh = 6. * z * (z - 1.) * field[i] / length +
            (1. - 4. * z + 3. * z2) * slope[i] +
            6. * z * (1. - z) * field[i + 1] / length +
            z * (3. * z - 2.) * slope[i + 1];
        v = h / b;
        dv = 0.5 * (dh / (b * b) - h / (b * b * b));
        return;
    }
    v = slope[last];
    dv = 0.;
}

__global__ void materialKernel(const DeviceElement *elements,
                               const DeviceMaterial *materials,
                               std::size_t elementCount,
                               const double *solution,
                               const double *bhFlux,
                               const double *bhField,
                               const double *bhSlope,
                               int iteration, bool warmStart,
                               double *mu1, double *mu2, double *dv)
{
    const std::size_t index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= elementCount) return;
    const DeviceElement &element = elements[index];
    const DeviceMaterial &material = materials[element.material];
    const double fill = material.fill;

    if (iteration == 0 && !warmStart) {
        if (material.laminationType == 0) {
            mu1[index] = material.muX * fill + (1. - fill);
            mu2[index] = material.muY * fill + (1. - fill);
        } else if (material.laminationType == 1) {
            const double initial = material.muX;
            mu1[index] = initial * fill + (1. - fill);
            mu2[index] = initial / (fill + initial * (1. - fill));
        } else if (material.laminationType == 2) {
            const double initial = material.muY;
            mu2[index] = initial * fill + (1. - fill);
            mu1[index] = initial / (fill + initial * (1. - fill));
        } else {
            mu1[index] = 1.;
            mu2[index] = 1.;
        }
        dv[index] = 0.;
        return;
    }

    if (material.bhCount <= 0) {
        if (material.laminationType == 0) {
            mu1[index] = material.muX * fill + (1. - fill);
            mu2[index] = material.muY * fill + (1. - fill);
        } else if (material.laminationType == 1) {
            const double initial = material.muX;
            mu1[index] = initial * fill + (1. - fill);
            mu2[index] = initial / (fill + initial * (1. - fill));
        } else if (material.laminationType == 2) {
            const double initial = material.muY;
            mu2[index] = initial * fill + (1. - fill);
            mu1[index] = initial / (fill + initial * (1. - fill));
        } else {
            mu1[index] = 1.;
            mu2[index] = 1.;
        }
        dv[index] = 0.;
        return;
    }

    double b1 = 0.;
    double b2 = 0.;
    for (int local = 0; local < 3; ++local) {
        const double a = solution[element.nodes[local]];
        if (material.laminationType == 1) {
            b1 += a * element.q[local];
            b2 += a * element.p[local] / fill;
        } else if (material.laminationType == 2) {
            b1 += a * element.q[local] / fill;
            b2 += a * element.p[local];
        } else {
            b1 += a * element.q[local];
            b2 += a * element.p[local];
        }
    }
    const double B = kEquationScale * sqrt(b1 * b1 + b2 * b2) /
        (0.02 * element.area);
    double reluctivity = 0.;
    double differential = 0.;
    bhProperties(material, B, bhFlux, bhField, bhSlope,
                 reluctivity, differential);
    const double relativeMu = 1. / (kMu0 * reluctivity);
    if (material.laminationType == 1) {
        mu1[index] = relativeMu * fill;
        mu2[index] = relativeMu / (fill + relativeMu * (1. - fill));
    } else if (material.laminationType == 2) {
        mu2[index] = relativeMu * fill;
        mu1[index] = relativeMu / (fill + relativeMu * (1. - fill));
    } else {
        mu1[index] = relativeMu;
        mu2[index] = relativeMu;
    }
    dv[index] = differential;
}

__global__ void bhSampleKernel(DeviceMaterial material, const double *flux,
                               const double *field, const double *slope,
                               const double *samples, std::size_t count,
                               CudaPlanarBhResult *results)
{
    const std::size_t index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= count) return;
    bhProperties(material, samples[index], flux, field, slope,
                 results[index].reluctivity,
                 results[index].differentialReluctivity);
}

__global__ void elementKernel(const DeviceElement *elements,
                              const DeviceMaterial *materials,
                              std::size_t elementCount,
                              const double *solution,
                              const double *circuitSource,
                              const std::int32_t *circuitCase,
                              std::size_t circuitCount,
                              const double *mu1, const double *mu2,
                              const double *dv,
                              int iteration, bool warmStart,
                              double *localMatrix, double *localRhs)
{
    const std::size_t index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= elementCount) return;
    const DeviceElement &element = elements[index];
    const DeviceMaterial &material = materials[element.material];
    double nonlinear[9] = {};

    if ((iteration > 0 || warmStart) && material.bhCount > 0) {
        double v[3] = {};
        double u[3] = {};
        if (material.laminationType == 0 && material.muX == material.muY) {
            for (int row = 0; row < 3; ++row)
                for (int column = 0; column < 3; ++column)
                    v[row] += (element.mx[3 * row + column] +
                               element.my[3 * row + column]) *
                              solution[element.nodes[column]];
            const double coefficient = -200. * kEquationScale *
                kEquationScale * kEquationScale * dv[index] / element.area;
            for (int row = 0; row < 3; ++row)
                for (int column = 0; column < 3; ++column)
                    nonlinear[3 * row + column] =
                        coefficient * v[row] * v[column];
        } else if (material.laminationType == 1) {
            const double fill = material.fill;
            for (int row = 0; row < 3; ++row)
                for (int column = 0; column < 3; ++column) {
                    v[row] += (element.my[3 * row + column] / fill +
                               element.mx[3 * row + column]) *
                              solution[element.nodes[column]];
                    u[row] += (element.my[3 * row + column] / fill +
                               fill * element.mx[3 * row + column]) *
                              solution[element.nodes[column]];
                }
            const double coefficient = -100. * kEquationScale *
                kEquationScale * kEquationScale * dv[index] / element.area;
            for (int row = 0; row < 3; ++row)
                for (int column = 0; column < 3; ++column)
                    nonlinear[3 * row + column] = coefficient *
                        (v[row] * u[column] + v[column] * u[row]);
        } else if (material.laminationType == 2) {
            const double fill = material.fill;
            for (int row = 0; row < 3; ++row)
                for (int column = 0; column < 3; ++column) {
                    v[row] += (element.mx[3 * row + column] / fill +
                               element.my[3 * row + column]) *
                              solution[element.nodes[column]];
                    u[row] += (element.mx[3 * row + column] / fill +
                               fill * element.my[3 * row + column]) *
                              solution[element.nodes[column]];
                }
            const double coefficient = -100. * kEquationScale *
                kEquationScale * kEquationScale * dv[index] / element.area;
            for (int row = 0; row < 3; ++row)
                for (int column = 0; column < 3; ++column)
                    nonlinear[3 * row + column] = coefficient *
                        (v[row] * u[column] + v[column] * u[row]);
        }
    }

    double source = 0.;
    if (element.circuit >= 0 &&
        static_cast<std::size_t>(element.circuit) < circuitCount) {
        source = circuitSource[element.circuit];
        if (circuitCase[element.circuit] == 0)
            source *= material.conductivity;
    }
    const double distributedRhs = -source * element.area / 3.;
    double rhs[3];
    double matrix[9];
    for (int row = 0; row < 3; ++row) {
        rhs[row] = element.fixedRhs[row] + distributedRhs;
        for (int column = 0; column < 3; ++column) {
            const int slot = 3 * row + column;
            matrix[slot] = element.fixedMatrix[slot] +
                element.mx[slot] / mu2[index] +
                element.my[slot] / mu1[index] + nonlinear[slot];
            rhs[row] += nonlinear[slot] * solution[element.nodes[column]];
        }
    }
    int upper = 0;
    for (int row = 0; row < 3; ++row)
        for (int column = row; column < 3; ++column)
            localMatrix[index * 6 + upper++] = -matrix[3 * row + column];
    for (int row = 0; row < 3; ++row)
        localRhs[index * 3 + row] = -rhs[row];
}

__global__ void atomicScatterKernel(const std::int32_t *matrixDestinations,
                                    const std::int32_t *rhsDestinations,
                                    const double *localMatrix,
                                    const double *localRhs,
                                    std::size_t elementCount,
                                    double *values, double *rhs)
{
    const std::size_t index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= elementCount) return;
    for (int local = 0; local < 6; ++local)
        atomicAdd(values + matrixDestinations[index * 6 + local],
                  localMatrix[index * 6 + local]);
    for (int local = 0; local < 3; ++local)
        atomicAdd(rhs + rhsDestinations[index * 3 + local],
                  localRhs[index * 3 + local]);
}

__global__ void deterministicMatrixScatterKernel(
    const std::int32_t *offsets, const std::int32_t *references,
    const double *localMatrix, std::size_t entryCount, double *values)
{
    const std::size_t entry = blockIdx.x * blockDim.x + threadIdx.x;
    if (entry >= entryCount) return;
    double sum = 0.;
    for (std::int32_t i = offsets[entry]; i < offsets[entry + 1]; ++i)
        sum += localMatrix[references[i]];
    values[entry] = sum;
}

__global__ void deterministicRhsScatterKernel(
    const std::int32_t *offsets, const std::int32_t *references,
    const double *localRhs, const double *explicitRhs,
    std::size_t nodeCount, double *rhs)
{
    const std::size_t node = blockIdx.x * blockDim.x + threadIdx.x;
    if (node >= nodeCount) return;
    double sum = explicitRhs[node];
    for (std::int32_t i = offsets[node]; i < offsets[node + 1]; ++i)
        sum += localRhs[references[i]];
    rhs[node] = sum;
}

__global__ void initializeRhsKernel(double *rhs, const double *explicitRhs,
                                    std::size_t nodeCount)
{
    const std::size_t node = blockIdx.x * blockDim.x + threadIdx.x;
    if (node < nodeCount) rhs[node] = explicitRhs[node];
}

__global__ void addContributionsKernel(const DeviceAdd *adds,
                                       std::size_t count, double *values,
                                       bool atomic)
{
    const std::size_t index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= count) return;
    if (atomic) atomicAdd(values + adds[index].destination, adds[index].value);
    else values[adds[index].destination] += adds[index].value;
}

__device__ std::int32_t findEntry(const std::int32_t *rows,
                                  const std::int32_t *columns,
                                  std::int32_t first, std::int32_t second)
{
    const std::int32_t row = max(first, second);
    const std::int32_t column = min(first, second);
    std::int32_t low = rows[row];
    std::int32_t high = rows[row + 1];
    while (low < high) {
        const std::int32_t middle = low + (high - low) / 2;
        if (columns[middle] < column) low = middle + 1;
        else high = middle;
    }
    return low < rows[row + 1] && columns[low] == column ? low : -1;
}

__global__ void dirichletKernel(const std::int32_t *rows,
                                const std::int32_t *columns,
                                std::size_t entryCount, std::int32_t node,
                                double value, double *values, double *rhs)
{
    const std::size_t entry = blockIdx.x * blockDim.x + threadIdx.x;
    if (entry < entryCount) {
        // The caller launches one thread per row instead; entry is a row here.
        const std::int32_t row = static_cast<std::int32_t>(entry);
        for (std::int32_t position = rows[row]; position < rows[row + 1]; ++position) {
            const std::int32_t column = columns[position];
            if (row == node && column != node) {
                rhs[column] -= values[position] * value;
                values[position] = 0.;
            } else if (column == node && row != node) {
                rhs[row] -= values[position] * value;
                values[position] = 0.;
            }
        }
    }
    if (entry == 0) {
        const std::int32_t diagonal = findEntry(rows, columns, node, node);
        rhs[node] = diagonal >= 0 ? values[diagonal] * value : 0.;
    }
}

__global__ void periodicKernel(const std::int32_t *rows,
                               const std::int32_t *columns,
                               std::size_t nodeCount, std::int32_t first,
                               std::int32_t second, bool anti,
                               double *values, double *rhs)
{
    const std::size_t k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < nodeCount && k != static_cast<std::size_t>(first) &&
        k != static_cast<std::size_t>(second)) {
        const std::int32_t left = findEntry(rows, columns,
                                            static_cast<std::int32_t>(k), first);
        const std::int32_t right = findEntry(rows, columns,
                                             static_cast<std::int32_t>(k), second);
        if (left >= 0 && right >= 0) {
            const double combined = anti
                ? (values[left] - values[right]) / 2.
                : (values[left] + values[right]) / 2.;
            values[left] = combined;
            values[right] = anti ? -combined : combined;
        }
    }
    if (k == 0) {
        const std::int32_t left = findEntry(rows, columns, first, first);
        const std::int32_t right = findEntry(rows, columns, second, second);
        const double diagonal = (values[left] + values[right]) / 2.;
        values[left] = diagonal;
        values[right] = diagonal;
        const double combined = anti ? (rhs[first] - rhs[second]) / 2.
                                     : (rhs[first] + rhs[second]) / 2.;
        rhs[first] = combined;
        rhs[second] = anti ? -combined : combined;
    }
}

__global__ void symmetricProductKernel(const std::int32_t *rows,
                                       const std::int32_t *columns,
                                       const double *values,
                                       const double *solution,
                                       std::size_t nodeCount,
                                       double *product)
{
    const std::size_t row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row >= nodeCount) return;
    for (std::int32_t entry = rows[row]; entry < rows[row + 1]; ++entry) {
        const std::int32_t column = columns[entry];
        const double value = values[entry];
        atomicAdd(product + row, value * solution[column]);
        if (column != static_cast<std::int32_t>(row))
            atomicAdd(product + column, value * solution[row]);
    }
}

__global__ void residualNormKernel(const double *product, const double *rhs,
                                   std::size_t nodeCount, double *sums)
{
    const std::size_t node = blockIdx.x * blockDim.x + threadIdx.x;
    if (node >= nodeCount) return;
    const double residual = product[node] - rhs[node];
    atomicAdd(sums, residual * residual);
    atomicAdd(sums + 1, rhs[node] * rhs[node]);
}

constexpr int kThreads = 256;
std::size_t blocks(std::size_t count)
{
    return (count + kThreads - 1) / kThreads;
}

} // namespace

class CudaPlanarAssembly::Implementation
{
public:
    Implementation(const PlanarAssemblyPlan &plan,
                   const std::vector<std::int32_t> &rows,
                   const std::vector<std::int32_t> &columns,
                   PlanarAssemblyBackend requested)
        : hostRows(rows), hostColumns(columns), mode(requested),
          nodeCount(static_cast<std::size_t>(plan.nodeCount)),
          elementCount(plan.elements.size()), entryCount(columns.size())
    {
        if (mode == PlanarAssemblyBackend::Host)
            throw std::invalid_argument("host mode cannot construct CUDA planar assembly");
        if (nodeCount + 1 != rows.size() || plan.explicitNodalRhs.size() != nodeCount)
            throw std::invalid_argument("CUDA planar plan dimensions disagree with CSR");

        std::vector<DeviceElement> elements(elementCount);
        std::vector<std::int32_t> matrixDestinations(elementCount * 6);
        std::vector<std::int32_t> rhsDestinations(elementCount * 3);
        std::vector<std::vector<std::int32_t>> matrixReferences(entryCount);
        std::vector<std::vector<std::int32_t>> rhsReferences(nodeCount);
        std::size_t circuitCount = 0;
        for (std::size_t i = 0; i < elementCount; ++i) {
            const auto &source = plan.elements[i];
            auto &target = elements[i];
            std::copy(source.nodes.begin(), source.nodes.end(), target.nodes);
            std::copy(source.p.begin(), source.p.end(), target.p);
            std::copy(source.q.begin(), source.q.end(), target.q);
            std::copy(source.mx.begin(), source.mx.end(), target.mx);
            std::copy(source.my.begin(), source.my.end(), target.my);
            std::copy(source.mxy.begin(), source.mxy.end(), target.mxy);
            std::copy(source.fixedMatrix.begin(), source.fixedMatrix.end(),
                      target.fixedMatrix);
            std::copy(source.fixedRhs.begin(), source.fixedRhs.end(), target.fixedRhs);
            target.area = source.area;
            target.material = source.material;
            target.circuit = source.circuit;
            if (source.circuit >= 0)
                circuitCount = std::max(circuitCount,
                    static_cast<std::size_t>(source.circuit) + 1);
            int upper = 0;
            for (int row = 0; row < 3; ++row) {
                rhsDestinations[i * 3 + row] = source.nodes[row];
                rhsReferences[static_cast<std::size_t>(source.nodes[row])]
                    .push_back(static_cast<std::int32_t>(i * 3 + row));
                for (int column = row; column < 3; ++column) {
                    const std::int32_t destination = findLowerEntry(
                        rows, columns, source.nodes[row], source.nodes[column]);
                    matrixDestinations[i * 6 + upper] = destination;
                    matrixReferences[static_cast<std::size_t>(destination)]
                        .push_back(static_cast<std::int32_t>(i * 6 + upper));
                    ++upper;
                }
            }
        }
        this->circuitCount = circuitCount;

        std::vector<DeviceMaterial> materials(plan.materials.size());
        for (std::size_t i = 0; i < materials.size(); ++i) {
            materials[i] = {plan.materials[i].muX, plan.materials[i].muY,
                            plan.materials[i].laminationFill,
                            plan.materials[i].conductivity,
                            plan.materials[i].laminationType,
                            plan.materials[i].bhOffset,
                            plan.materials[i].bhCount};
        }

        flatten(matrixReferences, matrixOffsetsHost, matrixReferencesHost);
        flatten(rhsReferences, rhsOffsetsHost, rhsReferencesHost);
        dRows.resize(rows.size());
        dColumns.resize(columns.size());
        dElements.resize(elements.size());
        dMaterials.resize(materials.size());
        dBhFlux.resize(plan.bhFluxDensity.size());
        dBhField.resize(plan.bhField.size());
        dBhSlope.resize(plan.bhSlope.size());
        dExplicitRhs.resize(plan.explicitNodalRhs.size());
        dMatrixDestinations.resize(matrixDestinations.size());
        dRhsDestinations.resize(rhsDestinations.size());
        dMatrixOffsets.resize(matrixOffsetsHost.size());
        dMatrixReferences.resize(matrixReferencesHost.size());
        dRhsOffsets.resize(rhsOffsetsHost.size());
        dRhsReferences.resize(rhsReferencesHost.size());
        dCircuitSource.resize(circuitCount);
        dCircuitCase.resize(circuitCount);
        dMu1.resize(elementCount);
        dMu2.resize(elementCount);
        dDv.resize(elementCount);
        dLocalMatrix.resize(elementCount * 6);
        dLocalRhs.resize(elementCount * 3);
        dResidualProduct.resize(nodeCount);
        dResidualSums.resize(2);

        copyHost(dRows, rows);
        copyHost(dColumns, columns);
        copyHost(dElements, elements);
        copyHost(dMaterials, materials);
        copyHost(dBhFlux, plan.bhFluxDensity);
        copyHost(dBhField, plan.bhField);
        copyHost(dBhSlope, plan.bhSlope);
        copyHost(dExplicitRhs, plan.explicitNodalRhs);
        copyHost(dMatrixDestinations, matrixDestinations);
        copyHost(dRhsDestinations, rhsDestinations);
        copyHost(dMatrixOffsets, matrixOffsetsHost);
        copyHost(dMatrixReferences, matrixReferencesHost);
        copyHost(dRhsOffsets, rhsOffsetsHost);
        copyHost(dRhsReferences, rhsReferencesHost);
    }

    template <typename T>
    static void copyHost(Buffer<T> &destination, const std::vector<T> &source)
    {
        if (!source.empty())
            XFEMM_CUDA_ASSEMBLY_CHECK(cudaMemcpy(
                destination.get(), source.data(), destination.bytes(),
                cudaMemcpyHostToDevice));
    }

    static void flatten(const std::vector<std::vector<std::int32_t>> &source,
                        std::vector<std::int32_t> &offsets,
                        std::vector<std::int32_t> &references)
    {
        offsets.resize(source.size() + 1);
        for (std::size_t i = 0; i < source.size(); ++i) {
            offsets[i] = static_cast<std::int32_t>(references.size());
            references.insert(references.end(), source[i].begin(), source[i].end());
        }
        offsets.back() = static_cast<std::int32_t>(references.size());
    }

    std::vector<DeviceAdd> mapAge(const std::vector<CudaPlanarMatrixAdd> &source,
                                  bool deterministic) const
    {
        if (deterministic) {
            std::map<std::int32_t, double> reduced;
            for (const auto &entry : source)
                reduced[findLowerEntry(hostRows, hostColumns,
                                       entry.row, entry.column)] += entry.value;
            std::vector<DeviceAdd> result;
            result.reserve(reduced.size());
            for (const auto &entry : reduced)
                result.push_back({entry.first, entry.second});
            return result;
        }
        std::vector<DeviceAdd> result;
        result.reserve(source.size());
        for (const auto &entry : source)
            result.push_back({findLowerEntry(hostRows, hostColumns,
                                             entry.row, entry.column), entry.value});
        return result;
    }

    std::vector<std::int32_t> hostRows;
    std::vector<std::int32_t> hostColumns;
    PlanarAssemblyBackend mode;
    std::size_t nodeCount = 0;
    std::size_t elementCount = 0;
    std::size_t entryCount = 0;
    std::size_t circuitCount = 0;
    std::vector<std::int32_t> matrixOffsetsHost;
    std::vector<std::int32_t> matrixReferencesHost;
    std::vector<std::int32_t> rhsOffsetsHost;
    std::vector<std::int32_t> rhsReferencesHost;
    Buffer<std::int32_t> dRows, dColumns;
    Buffer<DeviceElement> dElements;
    Buffer<DeviceMaterial> dMaterials;
    Buffer<double> dBhFlux, dBhField, dBhSlope, dExplicitRhs;
    Buffer<std::int32_t> dMatrixDestinations, dRhsDestinations;
    Buffer<std::int32_t> dMatrixOffsets, dMatrixReferences;
    Buffer<std::int32_t> dRhsOffsets, dRhsReferences;
    Buffer<double> dCircuitSource, dMu1, dMu2, dDv;
    Buffer<std::int32_t> dCircuitCase;
    Buffer<double> dLocalMatrix, dLocalRhs;
    Buffer<DeviceAdd> dAgeAdds;
    Buffer<double> dResidualProduct, dResidualSums;
    std::array<Event, 6> assemblyEvents;
    Event constraintBegin;
    Event constraintEnd;

    std::uint64_t deviceBytes() const
    {
        return dRows.bytes() + dColumns.bytes() + dElements.bytes() +
            dMaterials.bytes() + dBhFlux.bytes() + dBhField.bytes() +
            dBhSlope.bytes() + dExplicitRhs.bytes() +
            dMatrixDestinations.bytes() + dRhsDestinations.bytes() +
            dMatrixOffsets.bytes() + dMatrixReferences.bytes() +
            dRhsOffsets.bytes() + dRhsReferences.bytes() +
            dCircuitSource.bytes() + dCircuitCase.bytes() + dMu1.bytes() +
            dMu2.bytes() + dDv.bytes() + dLocalMatrix.bytes() +
            dLocalRhs.bytes() + dAgeAdds.bytes() + dResidualProduct.bytes() +
            dResidualSums.bytes();
    }
};

CudaPlanarAssembly::CudaPlanarAssembly(
    const PlanarAssemblyPlan &plan,
    const std::vector<std::int32_t> &lowerRows,
    const std::vector<std::int32_t> &lowerColumns,
    PlanarAssemblyBackend mode)
    : m_impl(std::make_unique<Implementation>(plan, lowerRows, lowerColumns, mode))
{}

CudaPlanarAssembly::~CudaPlanarAssembly() = default;

std::vector<CudaPlanarBhResult> evaluateCudaPlanarBh(
    const PlanarAssemblyMaterial &material,
    const std::vector<double> &bhFluxDensity,
    const std::vector<double> &bhField,
    const std::vector<double> &bhSlope,
    const std::vector<double> &fluxDensitySamples)
{
    if (bhFluxDensity.size() != bhField.size() ||
        bhFluxDensity.size() != bhSlope.size() ||
        material.bhOffset < 0 || material.bhCount < 0 ||
        static_cast<std::size_t>(material.bhOffset + material.bhCount) >
            bhFluxDensity.size())
        throw std::invalid_argument("CUDA B-H sample dimensions disagree");
    DeviceMaterial deviceMaterial{material.muX, material.muY,
        material.laminationFill, material.conductivity,
        material.laminationType, material.bhOffset, material.bhCount};
    Buffer<double> dFlux(bhFluxDensity.size());
    Buffer<double> dField(bhField.size());
    Buffer<double> dSlope(bhSlope.size());
    Buffer<double> dSamples(fluxDensitySamples.size());
    Buffer<CudaPlanarBhResult> dResults(fluxDensitySamples.size());
    if (!bhFluxDensity.empty()) {
        XFEMM_CUDA_ASSEMBLY_CHECK(cudaMemcpy(
            dFlux.get(), bhFluxDensity.data(), dFlux.bytes(),
            cudaMemcpyHostToDevice));
        XFEMM_CUDA_ASSEMBLY_CHECK(cudaMemcpy(
            dField.get(), bhField.data(), dField.bytes(),
            cudaMemcpyHostToDevice));
        XFEMM_CUDA_ASSEMBLY_CHECK(cudaMemcpy(
            dSlope.get(), bhSlope.data(), dSlope.bytes(),
            cudaMemcpyHostToDevice));
    }
    if (!fluxDensitySamples.empty()) {
        XFEMM_CUDA_ASSEMBLY_CHECK(cudaMemcpy(
            dSamples.get(), fluxDensitySamples.data(), dSamples.bytes(),
            cudaMemcpyHostToDevice));
        bhSampleKernel<<<blocks(fluxDensitySamples.size()), kThreads>>>(
            deviceMaterial, dFlux.get(), dField.get(), dSlope.get(),
            dSamples.get(), fluxDensitySamples.size(), dResults.get());
        XFEMM_CUDA_ASSEMBLY_CHECK(cudaGetLastError());
    }
    std::vector<CudaPlanarBhResult> results(fluxDensitySamples.size());
    if (!results.empty())
        XFEMM_CUDA_ASSEMBLY_CHECK(cudaMemcpy(
            results.data(), dResults.get(), dResults.bytes(),
            cudaMemcpyDeviceToHost));
    return results;
}

CudaPlanarAssemblyTimings CudaPlanarAssembly::assemble(
    void *streamPointer, double *deviceValues, double *deviceRhs,
    double *deviceSolution, const double *hostSolution,
    const PlanarAssemblyState &state,
    const std::vector<CudaPlanarMatrixAdd> &ageContributions)
{
    auto &p = *m_impl;
    auto stream = static_cast<cudaStream_t>(streamPointer);
    if (state.circuitSource.size() != p.circuitCount)
        throw std::invalid_argument("CUDA planar circuit-source count mismatch");
    if (state.circuitCase.size() != p.circuitCount)
        throw std::invalid_argument("CUDA planar circuit-case count mismatch");
    CudaPlanarAssemblyTimings timings;
    const bool atomic = p.mode == PlanarAssemblyBackend::CudaAtomic;
    std::vector<DeviceAdd> age = p.mapAge(ageContributions, !atomic);
    p.dAgeAdds.resize(age.size());

    XFEMM_CUDA_ASSEMBLY_CHECK(
        cudaEventRecord(p.assemblyEvents[0].get(), stream));
    XFEMM_CUDA_ASSEMBLY_CHECK(cudaMemcpyAsync(
        deviceSolution, hostSolution, p.nodeCount * sizeof(double),
        cudaMemcpyHostToDevice, stream));
    timings.transferBytes += p.nodeCount * sizeof(double);
    if (p.circuitCount) {
        XFEMM_CUDA_ASSEMBLY_CHECK(cudaMemcpyAsync(
            p.dCircuitSource.get(), state.circuitSource.data(),
            p.dCircuitSource.bytes(), cudaMemcpyHostToDevice, stream));
        timings.transferBytes += p.dCircuitSource.bytes();
        XFEMM_CUDA_ASSEMBLY_CHECK(cudaMemcpyAsync(
            p.dCircuitCase.get(), state.circuitCase.data(),
            p.dCircuitCase.bytes(), cudaMemcpyHostToDevice, stream));
        timings.transferBytes += p.dCircuitCase.bytes();
    }
    if (atomic) {
        XFEMM_CUDA_ASSEMBLY_CHECK(cudaMemsetAsync(
            deviceValues, 0, p.entryCount * sizeof(double), stream));
        initializeRhsKernel<<<blocks(p.nodeCount), kThreads, 0, stream>>>(
            deviceRhs, p.dExplicitRhs.get(), p.nodeCount);
        XFEMM_CUDA_ASSEMBLY_CHECK(cudaGetLastError());
    }
    XFEMM_CUDA_ASSEMBLY_CHECK(
        cudaEventRecord(p.assemblyEvents[1].get(), stream));

    materialKernel<<<blocks(p.elementCount), kThreads, 0, stream>>>(
        p.dElements.get(), p.dMaterials.get(), p.elementCount,
        deviceSolution, p.dBhFlux.get(), p.dBhField.get(), p.dBhSlope.get(),
        state.nonlinearIteration, state.warmStart,
        p.dMu1.get(), p.dMu2.get(), p.dDv.get());
    XFEMM_CUDA_ASSEMBLY_CHECK(cudaGetLastError());
    XFEMM_CUDA_ASSEMBLY_CHECK(
        cudaEventRecord(p.assemblyEvents[2].get(), stream));

    elementKernel<<<blocks(p.elementCount), kThreads, 0, stream>>>(
        p.dElements.get(), p.dMaterials.get(), p.elementCount,
        deviceSolution, p.dCircuitSource.get(), p.dCircuitCase.get(),
        p.circuitCount, p.dMu1.get(), p.dMu2.get(), p.dDv.get(),
        state.nonlinearIteration, state.warmStart,
        p.dLocalMatrix.get(), p.dLocalRhs.get());
    XFEMM_CUDA_ASSEMBLY_CHECK(cudaGetLastError());
    XFEMM_CUDA_ASSEMBLY_CHECK(
        cudaEventRecord(p.assemblyEvents[3].get(), stream));

    if (atomic) {
        atomicScatterKernel<<<blocks(p.elementCount), kThreads, 0, stream>>>(
            p.dMatrixDestinations.get(), p.dRhsDestinations.get(),
            p.dLocalMatrix.get(), p.dLocalRhs.get(), p.elementCount,
            deviceValues, deviceRhs);
        XFEMM_CUDA_ASSEMBLY_CHECK(cudaGetLastError());
    } else {
        deterministicMatrixScatterKernel<<<blocks(p.entryCount), kThreads, 0, stream>>>(
            p.dMatrixOffsets.get(), p.dMatrixReferences.get(),
            p.dLocalMatrix.get(), p.entryCount, deviceValues);
        XFEMM_CUDA_ASSEMBLY_CHECK(cudaGetLastError());
        deterministicRhsScatterKernel<<<blocks(p.nodeCount), kThreads, 0, stream>>>(
            p.dRhsOffsets.get(), p.dRhsReferences.get(), p.dLocalRhs.get(),
            p.dExplicitRhs.get(), p.nodeCount, deviceRhs);
        XFEMM_CUDA_ASSEMBLY_CHECK(cudaGetLastError());
    }
    XFEMM_CUDA_ASSEMBLY_CHECK(
        cudaEventRecord(p.assemblyEvents[4].get(), stream));

    if (!age.empty()) {
        XFEMM_CUDA_ASSEMBLY_CHECK(cudaMemcpyAsync(
            p.dAgeAdds.get(), age.data(), p.dAgeAdds.bytes(),
            cudaMemcpyHostToDevice, stream));
        timings.transferBytes += p.dAgeAdds.bytes();
        addContributionsKernel<<<blocks(age.size()), kThreads, 0, stream>>>(
            p.dAgeAdds.get(), age.size(), deviceValues, atomic);
        XFEMM_CUDA_ASSEMBLY_CHECK(cudaGetLastError());
    }
    XFEMM_CUDA_ASSEMBLY_CHECK(
        cudaEventRecord(p.assemblyEvents[5].get(), stream));
    XFEMM_CUDA_ASSEMBLY_CHECK(
        cudaEventSynchronize(p.assemblyEvents[5].get()));

    timings.clearMs = elapsedDevice(p.assemblyEvents[0], p.assemblyEvents[1]);
    timings.materialMs = elapsedDevice(p.assemblyEvents[1], p.assemblyEvents[2]);
    timings.elementMs = elapsedDevice(p.assemblyEvents[2], p.assemblyEvents[3]);
    timings.scatterMs = elapsedDevice(p.assemblyEvents[3], p.assemblyEvents[4]);
    timings.ageUploadMs = elapsedDevice(p.assemblyEvents[4], p.assemblyEvents[5]);
    return timings;
}

double CudaPlanarAssembly::applyConstraints(
    void *streamPointer, double *deviceValues, double *deviceRhs,
    const std::vector<CudaPlanarConstraint> &constraints)
{
    auto &p = *m_impl;
    auto stream = static_cast<cudaStream_t>(streamPointer);
    XFEMM_CUDA_ASSEMBLY_CHECK(cudaEventRecord(p.constraintBegin.get(), stream));
    for (const auto &constraint : constraints) {
        if (constraint.kind == CudaPlanarConstraintKind::Dirichlet) {
            dirichletKernel<<<blocks(p.nodeCount), kThreads, 0, stream>>>(
                p.dRows.get(), p.dColumns.get(), p.nodeCount,
                constraint.a, constraint.value, deviceValues, deviceRhs);
        } else {
            periodicKernel<<<blocks(p.nodeCount), kThreads, 0, stream>>>(
                p.dRows.get(), p.dColumns.get(), p.nodeCount,
                constraint.a, constraint.b,
                constraint.kind == CudaPlanarConstraintKind::Antiperiodic,
                deviceValues, deviceRhs);
        }
        XFEMM_CUDA_ASSEMBLY_CHECK(cudaGetLastError());
    }
    XFEMM_CUDA_ASSEMBLY_CHECK(cudaEventRecord(p.constraintEnd.get(), stream));
    XFEMM_CUDA_ASSEMBLY_CHECK(cudaEventSynchronize(p.constraintEnd.get()));
    return elapsedDevice(p.constraintBegin, p.constraintEnd);
}

double CudaPlanarAssembly::relativeResidual(
    void *streamPointer, const double *deviceValues, const double *deviceRhs,
    const double *deviceSolution)
{
    auto &p = *m_impl;
    auto stream = static_cast<cudaStream_t>(streamPointer);
    XFEMM_CUDA_ASSEMBLY_CHECK(cudaMemsetAsync(
        p.dResidualProduct.get(), 0, p.dResidualProduct.bytes(), stream));
    XFEMM_CUDA_ASSEMBLY_CHECK(cudaMemsetAsync(
        p.dResidualSums.get(), 0, p.dResidualSums.bytes(), stream));
    symmetricProductKernel<<<blocks(p.nodeCount), kThreads, 0, stream>>>(
        p.dRows.get(), p.dColumns.get(), deviceValues, deviceSolution,
        p.nodeCount, p.dResidualProduct.get());
    residualNormKernel<<<blocks(p.nodeCount), kThreads, 0, stream>>>(
        p.dResidualProduct.get(), deviceRhs, p.nodeCount,
        p.dResidualSums.get());
    double sums[2] = {};
    XFEMM_CUDA_ASSEMBLY_CHECK(cudaMemcpyAsync(
        sums, p.dResidualSums.get(), sizeof(sums), cudaMemcpyDeviceToHost,
        stream));
    XFEMM_CUDA_ASSEMBLY_CHECK(cudaStreamSynchronize(stream));
    return sums[1] == 0. ? std::sqrt(sums[0]) : std::sqrt(sums[0] / sums[1]);
}

void CudaPlanarAssembly::downloadMatrixAndRhs(
    void *streamPointer, const double *deviceValues, const double *deviceRhs,
    std::vector<double> &values, std::vector<double> &rhs) const
{
    const auto &p = *m_impl;
    auto stream = static_cast<cudaStream_t>(streamPointer);
    values.resize(p.entryCount);
    rhs.resize(p.nodeCount);
    XFEMM_CUDA_ASSEMBLY_CHECK(cudaMemcpyAsync(
        values.data(), deviceValues, values.size() * sizeof(double),
        cudaMemcpyDeviceToHost, stream));
    XFEMM_CUDA_ASSEMBLY_CHECK(cudaMemcpyAsync(
        rhs.data(), deviceRhs, rhs.size() * sizeof(double),
        cudaMemcpyDeviceToHost, stream));
    XFEMM_CUDA_ASSEMBLY_CHECK(cudaStreamSynchronize(stream));
}

std::uint64_t CudaPlanarAssembly::deviceBytes() const
{
    return m_impl->deviceBytes();
}

} // namespace femm
