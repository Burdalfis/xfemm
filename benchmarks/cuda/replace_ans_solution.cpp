#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::vector<double> readSolution(const std::string &path)
{
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("cannot open " + path);
    std::uint64_t count = 0;
    input.read(reinterpret_cast<char *>(&count), sizeof(count));
    std::vector<double> values(static_cast<std::size_t>(count));
    input.read(reinterpret_cast<char *>(values.data()),
               static_cast<std::streamsize>(values.size() * sizeof(double)));
    if (!input) throw std::runtime_error("truncated solution file " + path);
    return values;
}

std::string withoutCarriageReturn(std::string value)
{
    if (!value.empty() && value.back() == '\r') value.pop_back();
    return value;
}

} // namespace

int main(int argc, char **argv)
{
    if (argc != 4 && argc != 5) {
        std::cerr << "usage: xfemm-replace-ans-solution "
                     "TEMPLATE.ans SOLUTION.bin OUTPUT.ans [SCALE]\n";
        return 2;
    }
    try {
        const auto solution = readSolution(argv[2]);
        const double scale = argc == 5 ? std::stod(argv[4]) : 1.;
        std::ifstream input(argv[1]);
        if (!input) throw std::runtime_error(std::string("cannot open ") + argv[1]);
        std::ofstream output(argv[3], std::ios::trunc);
        if (!output) throw std::runtime_error(std::string("cannot create ") + argv[3]);

        std::string line;
        bool foundSolution = false;
        while (std::getline(input, line)) {
            output << line << '\n';
            if (withoutCarriageReturn(line) == "[Solution]") {
                foundSolution = true;
                break;
            }
        }
        if (!foundSolution || !std::getline(input, line))
            throw std::runtime_error("template has no [Solution] node section");
        output << line << '\n';
        const auto nodeCount = std::stoull(withoutCarriageReturn(line));
        if (nodeCount != solution.size())
            throw std::runtime_error("solution and template node counts differ");

        output << std::setprecision(17);
        for (std::size_t node = 0; node < solution.size(); ++node) {
            if (!std::getline(input, line))
                throw std::runtime_error("truncated template node section");
            std::istringstream fields(line);
            std::string x, y, ignoredPotential, boundary;
            if (!(fields >> x >> y >> ignoredPotential >> boundary))
                throw std::runtime_error("malformed template node record");
            output << x << '\t' << y << '\t' << solution[node] * scale
                   << '\t' << boundary;
            std::string extra;
            while (fields >> extra) output << '\t' << extra;
            output << '\n';
        }
        while (std::getline(input, line)) output << line << '\n';
        if (!output) throw std::runtime_error("failed while writing output .ans");
    } catch (const std::exception &error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
