#include "paw/ham_correction.hpp"

#include "paw/coulomb_correction.hpp"
#include "paw/kinetic_diff.hpp"

#include <limits>
#include <stdexcept>

namespace dft
{

namespace
{

std::size_t magnetic_channel_count(const std::vector<PAWState> &states)
{
    std::size_t count = 0;
    for (const PAWState &state : states)
    {
        count += static_cast<std::size_t>(2 * state.l + 1);
    }

    return count;
}

std::size_t magnetic_channel_count(const PAWSetup &setup)
{
    return magnetic_channel_count(setup.states());
}

std::size_t dense_matrix_row_count(const DenseMatrix &matrix)
{
    return static_cast<std::size_t>(matrix.shape()[0]);
}

std::size_t dense_matrix_column_count(const DenseMatrix &matrix)
{
    return static_cast<std::size_t>(matrix.shape()[1]);
}

int dense_matrix_extent(std::size_t value)
{
    if (value > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        throw std::overflow_error("DenseMatrix extent exceeds int range");
    }

    return static_cast<int>(value);
}

DenseMatrix make_dense_matrix(std::size_t row_count, std::size_t column_count)
{
    return DenseMatrix(typename DenseMatrix::ShapeType{dense_matrix_extent(row_count), dense_matrix_extent(column_count)});
}

bool is_square_dense_matrix(const DenseMatrix &matrix, std::size_t extent)
{
    return dense_matrix_row_count(matrix) == extent && dense_matrix_column_count(matrix) == extent;
}

DenseMatrix add_dense_matrices(const DenseMatrix &left, const DenseMatrix &right, const std::string &name)
{
    if (dense_matrix_row_count(left) != dense_matrix_row_count(right) ||
        dense_matrix_column_count(left) != dense_matrix_column_count(right))
    {
        throw std::runtime_error(name + ": matrix size mismatch");
    }

    DenseMatrix sum = left;
    for (std::size_t row_index = 0; row_index < dense_matrix_row_count(left); ++row_index)
    {
        for (std::size_t column_index = 0; column_index < dense_matrix_column_count(left); ++column_index)
        {
            const int row = dense_matrix_extent(row_index);
            const int column = dense_matrix_extent(column_index);
            sum[row, column] += right[row, column];
        }
    }

    return sum;
}

DenseMatrix expand_magnetic_matrix(const DenseMatrix &radial_matrix, const std::vector<PAWState> &states)
{
    const std::size_t radial_extent = dense_matrix_row_count(radial_matrix);
    if (radial_extent != states.size() || radial_extent != dense_matrix_column_count(radial_matrix))
    {
        throw std::runtime_error("Full radial correction matrix must match the number of states");
    }

    DenseMatrix expanded(
        typename DenseMatrix::ShapeType{static_cast<int>(magnetic_channel_count(states)),
                                        static_cast<int>(magnetic_channel_count(states))});
    std::vector<std::size_t> state_offset(states.size(), 0);

    std::size_t offset = 0;
    for (std::size_t state_index = 0; state_index < states.size(); ++state_index)
    {
        state_offset[state_index] = offset;
        offset += static_cast<std::size_t>(2 * states[state_index].l + 1);
    }

    for (std::size_t row_state_index = 0; row_state_index < states.size(); ++row_state_index)
    {
        for (int m = -states[row_state_index].l; m <= states[row_state_index].l; ++m)
        {
            const std::size_t row_index =
                state_offset[row_state_index] + static_cast<std::size_t>(m + states[row_state_index].l);
            for (std::size_t column_state_index = 0; column_state_index < states.size(); ++column_state_index)
            {
                if (std::abs(m) > states[column_state_index].l)
                {
                    continue;
                }

                const std::size_t column_index =
                    state_offset[column_state_index] + static_cast<std::size_t>(m + states[column_state_index].l);
                expanded[static_cast<int>(row_index), static_cast<int>(column_index)] =
                    radial_matrix[static_cast<int>(row_state_index), static_cast<int>(column_state_index)];
            }
        }
    }

    return expanded;
}

} // namespace

void HamCorrection::validate(const PAWSetup &setup) const
{
    const std::size_t full_extent = magnetic_channel_count(setup);
    if (!is_square_dense_matrix(m_kinetic_energy_differences, full_extent) ||
        !is_square_dense_matrix(m_static_coulomb_correction, full_extent) ||
        !is_square_dense_matrix(m_fixed_nonlocal_correction, full_extent))
    {
        throw std::runtime_error("Hamiltonian correction matrices must match the magnetic-channel count");
    }
}

HamCorrection build_ham_correction(const PAWSetup &setup)
{
    HamCorrection correction;
    const DenseMatrix radial_kinetic_energy_differences =
        build_full_kinetic_diff_matrix(setup.kinetic_difference_values(), setup.states().size());
    const DenseMatrix radial_static_coulomb_correction = build_two_index_coulomb_correction(setup);
    const DenseMatrix radial_fixed_nonlocal_correction = add_dense_matrices(radial_kinetic_energy_differences,
                                                                            radial_static_coulomb_correction,
                                                                "fixed_nonlocal_correction");
    correction.kinetic_energy_differences() = expand_magnetic_matrix(radial_kinetic_energy_differences, setup.states());
    correction.static_coulomb_correction() = expand_magnetic_matrix(radial_static_coulomb_correction, setup.states());
    correction.fixed_nonlocal_correction() = expand_magnetic_matrix(radial_fixed_nonlocal_correction, setup.states());
    correction.validate(setup);
    return correction;
}

} // namespace dft
