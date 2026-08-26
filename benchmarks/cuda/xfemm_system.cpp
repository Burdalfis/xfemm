#include "xfemm_system.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace xfemm_benchmark {
namespace {

template <typename T>
T readValue(std::ifstream &stream, const char *description)
{
    T value{};
    if (!stream.read(reinterpret_cast<char *>(&value), sizeof(value)))
        throw std::runtime_error(std::string("truncated system at ") + description);
    return value;
}

template <typename T>
void readArray(std::ifstream &stream, std::vector<T> &values,
               std::size_t count, const char *description)
{
    values.resize(count);
    if (count != 0 &&
        !stream.read(reinterpret_cast<char *>(values.data()),
                     static_cast<std::streamsize>(count * sizeof(T))))
        throw std::runtime_error(std::string("truncated system at ") + description);
}

double norm2(const std::vector<double> &values)
{
    long double sum = 0.;
    for (double value : values)
        sum += static_cast<long double>(value) * value;
    return std::sqrt(static_cast<double>(sum));
}

void hashBytes(std::uint64_t &hash, const void *data, std::size_t bytes)
{
    const auto *input = static_cast<const unsigned char *>(data);
    for (std::size_t i = 0; i < bytes; ++i) {
        hash ^= input[i];
        hash *= 1099511628211ull;
    }
}

} // namespace

LinearSystem readLinearSystem(const std::string &path)
{
    std::ifstream stream(path, std::ios::binary);
    if (!stream)
        throw std::runtime_error("cannot open system: " + path);

    char magic[8]{};
    if (!stream.read(magic, sizeof(magic)) ||
        std::memcmp(magic, "XFEMMLS\0", sizeof(magic)) != 0)
        throw std::runtime_error("not an XFEMM linear-system export: " + path);
    const std::uint32_t version = readValue<std::uint32_t>(stream, "version");
    const std::uint32_t endian = readValue<std::uint32_t>(stream, "endianness");
    const std::uint32_t scalarBytes = readValue<std::uint32_t>(stream, "scalar size");
    const std::uint32_t indexBytes = readValue<std::uint32_t>(stream, "index size");
    if (version != 1 || endian != 0x01020304u || scalarBytes != 8 || indexBytes != 4)
        throw std::runtime_error("unsupported XFEMM system format in: " + path);

    const std::uint64_t dimension64 = readValue<std::uint64_t>(stream, "dimension");
    const std::uint64_t nonzeros64 = readValue<std::uint64_t>(stream, "nonzero count");
    if (dimension64 > static_cast<std::uint64_t>(std::numeric_limits<int>::max()) ||
        nonzeros64 > static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
        throw std::runtime_error("system exceeds 32-bit CUDA sparse-library limits");
    const std::size_t dimension = static_cast<std::size_t>(dimension64);
    const std::size_t nonzeros = static_cast<std::size_t>(nonzeros64);

    LinearSystem result;
    result.sourcePath = path;
    result.solveIndex = readValue<std::uint64_t>(stream, "solve index");
    result.flags = readValue<std::uint32_t>(stream, "flags");
    result.cpuThreads = readValue<std::uint32_t>(stream, "CPU thread count");
    result.tolerance = readValue<double>(stream, "tolerance");
    result.ssorRelaxation = readValue<double>(stream, "SSOR relaxation");
    result.cpuIterations = readValue<std::int64_t>(stream, "CPU iterations");
    result.cpuPcgResidual = readValue<double>(stream, "CPU PCG residual");
    readArray(stream, result.rowOffsets, dimension + 1, "row offsets");
    readArray(stream, result.columnIndices, nonzeros, "column indices");
    readArray(stream, result.values, nonzeros, "matrix values");
    readArray(stream, result.rhs, dimension, "right-hand side");
    readArray(stream, result.initialSolution, dimension, "initial solution");
    readArray(stream, result.cpuSolution, dimension, "CPU solution");

    char trailing;
    if (stream.read(&trailing, 1))
        throw std::runtime_error("unexpected trailing data in: " + path);
    if (result.rowOffsets.empty() || result.rowOffsets.front() != 0 ||
        result.rowOffsets.back() != nonzeros)
        throw std::runtime_error("invalid CSR row offsets in: " + path);
    for (std::size_t i = 1; i < result.rowOffsets.size(); ++i)
        if (result.rowOffsets[i] < result.rowOffsets[i - 1])
            throw std::runtime_error("non-monotonic CSR row offsets in: " + path);
    for (std::int32_t column : result.columnIndices)
        if (column < 0 || static_cast<std::size_t>(column) >= dimension)
            throw std::runtime_error("out-of-range CSR column in: " + path);
    if ((result.flags & 2u) == 0)
        throw std::runtime_error("benchmark requires full symmetric CSR: " + path);
    return result;
}

std::uint64_t topologyHash(const LinearSystem &system)
{
    std::uint64_t hash = 1469598103934665603ull;
    hashBytes(hash, system.rowOffsets.data(),
              system.rowOffsets.size() * sizeof(system.rowOffsets[0]));
    hashBytes(hash, system.columnIndices.data(),
              system.columnIndices.size() * sizeof(system.columnIndices[0]));
    return hash;
}

std::uint64_t solutionHash(const std::vector<double> &solution)
{
    std::uint64_t hash = 1469598103934665603ull;
    hashBytes(hash, solution.data(), solution.size() * sizeof(solution[0]));
    return hash;
}

double relativeResidual(const LinearSystem &system,
                        const std::vector<double> &solution)
{
    if (solution.size() != system.dimension())
        throw std::invalid_argument("solution dimension mismatch");
    std::vector<double> residual(system.dimension());
    for (std::size_t row = 0; row < system.dimension(); ++row) {
        long double product = 0.;
        for (std::uint64_t j = system.rowOffsets[row];
             j < system.rowOffsets[row + 1]; ++j)
            product += static_cast<long double>(system.values[j]) *
                       solution[static_cast<std::size_t>(system.columnIndices[j])];
        residual[row] = system.rhs[row] - static_cast<double>(product);
    }
    const double denominator = norm2(system.rhs);
    return denominator == 0. ? norm2(residual) : norm2(residual) / denominator;
}

double relativeSolutionError(const std::vector<double> &candidate,
                             const std::vector<double> &reference)
{
    if (candidate.size() != reference.size())
        throw std::invalid_argument("solution dimension mismatch");
    std::vector<double> difference(candidate.size());
    std::transform(candidate.begin(), candidate.end(), reference.begin(),
                   difference.begin(), [](double a, double b) { return a - b; });
    const double denominator = norm2(reference);
    return denominator == 0. ? norm2(difference) : norm2(difference) / denominator;
}

} // namespace xfemm_benchmark
