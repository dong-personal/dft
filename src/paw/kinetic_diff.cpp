#include "paw/kinetic_diff.hpp"

#include <limits>
#include <stdexcept>

namespace dft
{

namespace
{

int kinetic_diff_extent(std::size_t value)
{
    if (value > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        throw std::overflow_error("Kinetic-difference matrix extent exceeds int range");
    }

    return static_cast<int>(value);
}

} // namespace

DenseMatrix build_full_kinetic_diff_matrix(const std::vector<double> &flat_values, std::size_t state_count)
{
    if (flat_values.size() != state_count * state_count)
    {
        throw std::runtime_error("kinetic_energy_differences: expected " + std::to_string(state_count * state_count) +
                                 " entries, got " + std::to_string(flat_values.size()));
    }

    DenseMatrix matrix(
        typename DenseMatrix::ShapeType{kinetic_diff_extent(state_count), kinetic_diff_extent(state_count)});
    for (std::size_t row_index = 0; row_index < state_count; ++row_index)
    {
        for (std::size_t column_index = 0; column_index < state_count; ++column_index)
        {
            matrix[kinetic_diff_extent(row_index), kinetic_diff_extent(column_index)] =
                flat_values[row_index * state_count + column_index];
        }
    }

    return matrix;
}

} // namespace dft
