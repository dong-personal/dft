#include "structure.h"

#include "constants.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace dft
{

namespace
{

constexpr int kSpatialDimension = 3;

std::string read_required_line(std::ifstream &input_stream, const std::string &error_context)
{
    std::string line;
    if (!std::getline(input_stream, line))
    {
        throw std::runtime_error("Failed to read POSCAR line for " + error_context);
    }
    return line;
}

std::vector<std::string> parse_element_symbols(const std::string &line)
{
    std::istringstream line_stream(line);
    std::vector<std::string> element_symbols;
    std::string element_symbol;
    while (line_stream >> element_symbol)
    {
        element_symbols.push_back(element_symbol);
    }
    return element_symbols;
}

std::vector<int> parse_element_counts(const std::string &line)
{
    std::istringstream line_stream(line);
    std::vector<int> element_counts;
    int element_count = 0;
    while (line_stream >> element_count)
    {
        element_counts.push_back(element_count);
    }
    return element_counts;
}

bool is_fractional_coordinate_line(const std::string &line)
{
    if (line.empty())
    {
        throw std::runtime_error("POSCAR coordinate type line is empty");
    }
    return line[0] == 'D' || line[0] == 'd';
}

Structure::LatticeVectors read_lattice(std::ifstream &input_stream, double scale)
{
    Structure::LatticeVectors lattice(typename Structure::LatticeVectors::ShapeType{3, 3});
    for (int row = 0; row < kSpatialDimension; ++row)
    {
        std::istringstream line_stream(read_required_line(input_stream, "lattice vector"));
        for (int col = 0; col < kSpatialDimension; ++col)
        {
            line_stream >> lattice[row, col];
            lattice[row, col] *= scale * ANG_TO_BOHR;
        }
    }
    return lattice;
}

Atom::AtomicPosition read_atomic_position(std::ifstream &input_stream, bool is_fractional_coordinate)
{
    Atom::AtomicPosition atomic_position(typename Atom::AtomicPosition::ShapeType{3});
    std::istringstream line_stream(read_required_line(input_stream, "atomic position"));
    for (int dim = 0; dim < kSpatialDimension; ++dim)
    {
        line_stream >> atomic_position[dim];
        if (!is_fractional_coordinate)
        {
            atomic_position[dim] *= ANG_TO_BOHR;
        }
    }
    return atomic_position;
}

} // namespace

Structure read_poscar(const std::string &filename)
{
    std::ifstream input_stream(filename);
    if (!input_stream)
    {
        throw std::runtime_error("Failed to open POSCAR: " + filename);
    }

    read_required_line(input_stream, "comment");
    const double scale = std::stod(read_required_line(input_stream, "scale"));
    const Structure::LatticeVectors lattice = read_lattice(input_stream, scale);

    const std::vector<std::string> element_symbols =
        parse_element_symbols(read_required_line(input_stream, "element symbols"));
    const std::vector<int> element_counts =
        parse_element_counts(read_required_line(input_stream, "element counts"));

    if (element_symbols.size() != element_counts.size())
    {
        throw std::runtime_error("POSCAR: element/count size mismatch");
    }

    const bool is_fractional_coordinate =
        is_fractional_coordinate_line(read_required_line(input_stream, "coordinate type"));

    Structure structure(
        lattice, is_fractional_coordinate ? Structure::CoordType::Fractional : Structure::CoordType::Cartesian);

    for (std::size_t element_index = 0; element_index < element_symbols.size(); ++element_index)
    {
        Atom atom(element_symbols[element_index], static_cast<std::size_t>(element_counts[element_index]));
        for (int atom_index = 0; atom_index < element_counts[element_index]; ++atom_index)
        {
            atom.set_position(static_cast<std::size_t>(atom_index),
                              read_atomic_position(input_stream, is_fractional_coordinate));
        }
        structure.add_atom(atom);
    }

    return structure;
}

} // namespace dft
