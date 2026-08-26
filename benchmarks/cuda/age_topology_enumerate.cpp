#include "xfemm_system.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace xb = xfemm_benchmark;

namespace {

std::vector<std::uint64_t> edges(const xb::LinearSystem &system)
{
    std::vector<std::uint64_t> result;
    for (std::size_t row = 0; row < system.dimension(); ++row)
        for (std::uint64_t j = system.rowOffsets[row];
             j < system.rowOffsets[row + 1]; ++j) {
            const auto column = static_cast<std::uint32_t>(system.columnIndices[j]);
            if (row < column)
                result.push_back((static_cast<std::uint64_t>(row) << 32) | column);
        }
    std::sort(result.begin(), result.end());
    return result;
}

std::vector<signed char> inferRingColors(const xb::LinearSystem &a,
                                         const xb::LinearSystem &b)
{
    const auto aEdges = edges(a);
    const auto bEdges = edges(b);
    std::vector<std::uint64_t> changed;
    std::set_symmetric_difference(aEdges.begin(), aEdges.end(), bEdges.begin(),
                                  bEdges.end(), std::back_inserter(changed));
    std::vector<std::vector<std::int32_t>> adjacency(a.dimension());
    for (const auto edge : changed) {
        const auto x = static_cast<std::int32_t>(edge >> 32);
        const auto y = static_cast<std::int32_t>(edge & 0xffffffffu);
        adjacency[x].push_back(y);
        adjacency[y].push_back(x);
    }
    std::vector<signed char> color(a.dimension(), -1);
    for (std::size_t start = 0; start < adjacency.size(); ++start) {
        if (adjacency[start].empty() || color[start] >= 0) continue;
        color[start] = 0;
        std::queue<std::int32_t> queue;
        queue.push(static_cast<std::int32_t>(start));
        while (!queue.empty()) {
            const auto node = queue.front();
            queue.pop();
            for (const auto neighbor : adjacency[node]) {
                if (color[neighbor] < 0) {
                    color[neighbor] = 1 - color[node];
                    queue.push(neighbor);
                } else if (color[neighbor] == color[node]) {
                    throw std::runtime_error("changed graph is not bipartite");
                }
            }
        }
    }
    return color;
}

std::vector<double> readNodeAngles(const std::string &path, std::size_t expected)
{
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open " + path);
    std::string line;
    while (std::getline(input, line) && line != "[Solution]" &&
           line != "[Solution]\r") {}
    if (!input || !std::getline(input, line))
        throw std::runtime_error(".ans has no [Solution] section");
    if (std::stoull(line) != expected)
        throw std::runtime_error(".ans and matrix dimensions differ");
    std::vector<double> angles(expected);
    const double twoPi = 2. * std::acos(-1.);
    for (std::size_t node = 0; node < expected; ++node) {
        if (!std::getline(input, line)) throw std::runtime_error("truncated .ans");
        std::istringstream fields(line);
        double x = 0., y = 0.;
        if (!(fields >> x >> y)) throw std::runtime_error("malformed .ans node");
        angles[node] = std::atan2(y, x);
        if (angles[node] < 0.) angles[node] += twoPi;
    }
    return angles;
}

std::uint64_t hashEdges(const std::vector<std::uint64_t> &edges)
{
    std::uint64_t hash = 1469598103934665603ull;
    for (const auto edge : edges) {
        for (unsigned byte = 0; byte < sizeof(edge); ++byte) {
            hash ^= static_cast<unsigned char>(edge >> (8 * byte));
            hash *= 1099511628211ull;
        }
    }
    return hash;
}

std::vector<std::uint64_t> crossEdges(const xb::LinearSystem &system,
                                      const std::vector<signed char> &color)
{
    auto all = edges(system);
    all.erase(std::remove_if(all.begin(), all.end(), [&](std::uint64_t edge) {
        const auto a = static_cast<std::uint32_t>(edge >> 32);
        const auto b = static_cast<std::uint32_t>(edge & 0xffffffffu);
        return color[a] < 0 || color[b] < 0 || color[a] == color[b];
    }), all.end());
    return all;
}

} // namespace

int main(int argc, char **argv)
{
    if (argc != 4) {
        std::cerr << "usage: xfemm-age-topology-enumerate "
                     "REFERENCE_SYSTEM SHIFTED_SYSTEM REFERENCE.ans\n";
        return 2;
    }
    try {
        const auto reference = xb::readLinearSystem(argv[1]);
        const auto shifted = xb::readLinearSystem(argv[2]);
        const auto color = inferRingColors(reference, shifted);
        const auto angles = readNodeAngles(argv[3], reference.dimension());
        std::vector<std::int32_t> first, second;
        for (std::size_t node = 0; node < color.size(); ++node) {
            if (color[node] == 0) first.push_back(static_cast<std::int32_t>(node));
            if (color[node] == 1) second.push_back(static_cast<std::int32_t>(node));
        }
        const auto byAngle = [&](std::int32_t a, std::int32_t b) {
            return angles[a] < angles[b];
        };
        std::sort(first.begin(), first.end(), byAngle);
        std::sort(second.begin(), second.end(), byAngle);
        if (first.size() != second.size() || first.empty())
            throw std::runtime_error("AGE rings do not have equal node counts");

        std::vector<std::size_t> secondPosition(reference.dimension());
        for (std::size_t i = 0; i < second.size(); ++i)
            secondPosition[second[i]] = i;
        const auto referenceCross = crossEdges(reference, color);
        const auto shiftedCross = crossEdges(shifted, color);
        const auto shiftedHash = hashEdges(shiftedCross);
        std::unordered_set<std::uint64_t> distinct;
        std::vector<std::size_t> shiftedMatches;
        for (std::size_t shift = 0; shift < second.size(); ++shift) {
            std::vector<std::uint64_t> transformed;
            transformed.reserve(referenceCross.size());
            for (auto edge : referenceCross) {
                std::uint32_t a = static_cast<std::uint32_t>(edge >> 32);
                std::uint32_t b = static_cast<std::uint32_t>(edge & 0xffffffffu);
                if (color[a] == 1) std::swap(a, b);
                b = static_cast<std::uint32_t>(
                    second[(secondPosition[b] + shift) % second.size()]);
                if (a > b) std::swap(a, b);
                transformed.push_back((static_cast<std::uint64_t>(a) << 32) | b);
            }
            std::sort(transformed.begin(), transformed.end());
            const auto hash = hashEdges(transformed);
            distinct.insert(hash);
            if (hash == shiftedHash && transformed == shiftedCross)
                shiftedMatches.push_back(shift);
        }
        std::cout << "ring_nodes," << first.size() << '\n'
                  << "reference_cross_edges," << referenceCross.size() << '\n'
                  << "shifted_cross_edges," << shiftedCross.size() << '\n'
                  << "distinct_cyclic_patterns," << distinct.size() << '\n'
                  << "mechanical_pattern_pitch_deg,"
                  << 360. / static_cast<double>(distinct.size()) << '\n'
                  << "shifted_pattern_matches";
        for (const auto match : shiftedMatches) std::cout << ',' << match;
        std::cout << '\n';
    } catch (const std::exception &error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
