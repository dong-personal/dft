#include "paw/grid_map.hpp"

#include <limits>
#include <stdexcept>

namespace dft
{

namespace
{

int grid_map_extent(std::size_t value)
{
    if (value > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        throw std::overflow_error("AtomGridMap extent exceeds int range");
    }

    return static_cast<int>(value);
}

DenseMatrix make_dense_matrix(std::size_t row_count, std::size_t column_count)
{
    return DenseMatrix(typename DenseMatrix::ShapeType{grid_map_extent(row_count), grid_map_extent(column_count)});
}

std::size_t row_count(const DenseMatrix &matrix)
{
    return static_cast<std::size_t>(matrix.shape()[0]);
}

std::size_t column_count(const DenseMatrix &matrix)
{
    return static_cast<std::size_t>(matrix.shape()[1]);
}

std::size_t projector_channel_count(const PAWSetup &setup)
{
    std::size_t channel_count = 0;
    for (const PAWState &state : setup.states())
    {
        channel_count += static_cast<std::size_t>(2 * state.l + 1);
    }

    return channel_count;
}

void validate_atom_gll_ham_correction_inputs(const DenseMatrix &projector_basis_overlap,
                                             const DenseMatrix &fixed_nonlocal, std::size_t dof_count)
{
    const std::size_t channel_count = row_count(projector_basis_overlap);

    if (column_count(projector_basis_overlap) != dof_count)
    {
        throw std::runtime_error("Projector-basis overlap column count must match the true-dof count");
    }

    if (row_count(fixed_nonlocal) != channel_count || column_count(fixed_nonlocal) != channel_count)
    {
        throw std::runtime_error("Fixed nonlocal correction size must match projector channel count");
    }
}

DenseMatrix assemble_atom_gll_ham_correction_matrix(const DenseMatrix &projector_basis_overlap,
                                                    const DenseMatrix &fixed_nonlocal, std::size_t dof_count)
{
    validate_atom_gll_ham_correction_inputs(projector_basis_overlap, fixed_nonlocal, dof_count);

    const std::size_t channel_count = row_count(projector_basis_overlap);
    DenseMatrix matrix = make_dense_matrix(dof_count, dof_count);
    for (std::size_t row_index = 0; row_index < dof_count; ++row_index)
    {
        for (std::size_t column_index = 0; column_index < dof_count; ++column_index)
        {
            double value = 0.0;
            for (std::size_t left_channel = 0; left_channel < channel_count; ++left_channel)
            {
                const double left_overlap =
                    projector_basis_overlap[grid_map_extent(left_channel), grid_map_extent(row_index)];
                if (left_overlap == 0.0)
                {
                    continue;
                }

                for (std::size_t right_channel = 0; right_channel < channel_count; ++right_channel)
                {
                    const double right_overlap =
                        projector_basis_overlap[grid_map_extent(right_channel), grid_map_extent(column_index)];
                    if (right_overlap == 0.0)
                    {
                        continue;
                    }

                    value += left_overlap *
                             fixed_nonlocal[grid_map_extent(left_channel), grid_map_extent(right_channel)] *
                             right_overlap;
                }
            }

            matrix[grid_map_extent(row_index), grid_map_extent(column_index)] = value;
        }
    }

    return matrix;
}

} // namespace

AtomGLLHamCorrectionMatrix::AtomGLLHamCorrectionMatrix(const DFTGLLHexSpace &space, const PAWSetup &setup,
                                                       const HamCorrection &ham_correction, const Atom &atom,
                                                       std::size_t position_index)
    : m_projector_basis_overlap(build_projector_basis_overlap_matrix(space, setup, atom, position_index)),
      m_fixed_nonlocal(ham_correction.fixed_nonlocal_correction()),
      m_matrix(assemble_atom_gll_ham_correction_matrix(m_projector_basis_overlap, m_fixed_nonlocal,
                                                       static_cast<std::size_t>(space.getDOF())))
{
    if (row_count(m_projector_basis_overlap) != projector_channel_count(setup))
    {
        throw std::runtime_error("Projector-basis overlap row count must match projector channel count");
    }
}

AtomGridMap build_atom_grid_map(const HamCorrection &ham_correction, const ProjectorTable &projector_table,
                                const std::vector<PeriodicNeighbor> &neighbors, const mfem::Vector &mass_diagonal)
{
    const DenseMatrix &fixed_nonlocal = ham_correction.fixed_nonlocal_correction();
    const std::size_t channel_count = projector_table.channel_count();

    if (row_count(projector_table.values) != channel_count)
    {
        throw std::runtime_error("Projector table row count must match projector channel metadata");
    }

    if (column_count(projector_table.values) != neighbors.size())
    {
        throw std::runtime_error("Projector table column count must match neighbor count");
    }

    if (row_count(fixed_nonlocal) != channel_count || column_count(fixed_nonlocal) != channel_count)
    {
        throw std::runtime_error("Fixed nonlocal correction size must match projector channel count");
    }

    AtomGridMap grid_map;
    grid_map.values = DenseMatrix(
        typename DenseMatrix::ShapeType{grid_map_extent(neighbors.size()), grid_map_extent(neighbors.size())});
    grid_map.point_indices.reserve(neighbors.size());

    for (const PeriodicNeighbor &neighbor : neighbors)
    {
        const int point_index = static_cast<int>(neighbor.index);
        if (point_index < 0 || point_index >= mass_diagonal.Size())
        {
            throw std::runtime_error("Neighbor point index exceeds mass-diagonal extent");
        }

        grid_map.point_indices.push_back(neighbor.index);
    }

    std::vector<double> weighted_projector(channel_count, 0.0);
    for (std::size_t source_index = 0; source_index < neighbors.size(); ++source_index)
    {
        const double source_mass = mass_diagonal(static_cast<int>(grid_map.point_indices[source_index]));

        for (std::size_t row_index = 0; row_index < channel_count; ++row_index)
        {
            double value = 0.0;
            for (std::size_t column_index = 0; column_index < channel_count; ++column_index)
            {
                value += fixed_nonlocal[grid_map_extent(row_index), grid_map_extent(column_index)] *
                         projector_table.values[grid_map_extent(column_index), grid_map_extent(source_index)];
            }

            weighted_projector[row_index] = source_mass * value;
        }

        for (std::size_t target_index = 0; target_index < neighbors.size(); ++target_index)
        {
            double value = 0.0;
            for (std::size_t projector_index = 0; projector_index < channel_count; ++projector_index)
            {
                value += projector_table.values[grid_map_extent(projector_index), grid_map_extent(target_index)] *
                         weighted_projector[projector_index];
            }

            grid_map.values[grid_map_extent(target_index), grid_map_extent(source_index)] = value;
        }
    }

    return grid_map;
}

DenseMatrix build_atom_gll_ham_correction_matrix(const DFTGLLHexSpace &space, const PAWSetup &setup,
                                                 const HamCorrection &ham_correction, const Atom &atom,
                                                 std::size_t position_index)
{
    const AtomGLLHamCorrectionMatrix atom_matrix(space, setup, ham_correction, atom, position_index);
    return atom_matrix.matrix();
}

DenseMatrix build_gll_ham_correction_matrix(const DFTGLLHexSpace &space, const PAWSetupRegistry &setup_registry)
{
    const Structure *structure = space.mesh_source().structure();
    if (structure == nullptr)
    {
        throw std::runtime_error(
            "DFTGLLHexSpace mesh source must have a structure before building GLL Hamiltonian correction");
    }

    DenseMatrix matrix =
        make_dense_matrix(static_cast<std::size_t>(space.getDOF()), static_cast<std::size_t>(space.getDOF()));

    for (const Atom &atom : structure->atoms())
    {
        const PAWSetup &setup = setup_registry.get(atom.symbol());
        const HamCorrection ham_correction = build_ham_correction(setup);

        for (std::size_t position_index = 0; position_index < atom.num_positions(); ++position_index)
        {
            const AtomGLLHamCorrectionMatrix atom_correction(space, setup, ham_correction, atom, position_index);
            const DenseMatrix &atom_matrix = atom_correction.matrix();

            for (std::size_t row_index = 0; row_index < row_count(matrix); ++row_index)
            {
                for (std::size_t column_index = 0; column_index < column_count(matrix); ++column_index)
                {
                    matrix[grid_map_extent(row_index), grid_map_extent(column_index)] +=
                        atom_matrix[grid_map_extent(row_index), grid_map_extent(column_index)];
                }
            }
        }
    }

    return matrix;
}

} // namespace dft
