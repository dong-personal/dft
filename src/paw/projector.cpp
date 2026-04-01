#include "paw/projector.hpp"

#include <array>
#include <cmath>
#include <limits>

namespace dft
{

namespace
{

using Point3 = Atom::AtomicPosition;

int projector_extent(std::size_t value)
{
    if (value > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        throw std::overflow_error("ProjectorValueTable extent exceeds int range");
    }

    return static_cast<int>(value);
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

Point3 make_point3(double x, double y, double z)
{
    Point3 point(typename Point3::ShapeType{3});
    point[0] = x;
    point[1] = y;
    point[2] = z;
    return point;
}

Point3 to_cartesian_point(const Atom::AtomicPosition &point)
{
    return make_point3(point[0], point[1], point[2]);
}

Point3 to_basis_point(const mfem::Vector &point)
{
    return make_point3(point[0], point[1], point[2]);
}

Point3 to_shifted_point(const Point3 &point, const Structure::LatticeVectors &lattice,
                       const std::array<int, 3> &periodic_image)
{
    return make_point3(point[0] + periodic_image[0] * lattice[0, 0] + periodic_image[1] * lattice[1, 0] +
                           periodic_image[2] * lattice[2, 0],
                       point[1] + periodic_image[0] * lattice[0, 1] + periodic_image[1] * lattice[1, 1] +
                           periodic_image[2] * lattice[2, 1],
                       point[2] + periodic_image[0] * lattice[0, 2] + periodic_image[1] * lattice[1, 2] +
                           periodic_image[2] * lattice[2, 2]);
}

Point3 subtract_vectors(const Point3 &left, const Point3 &right)
{
    return make_point3(left[0] - right[0], left[1] - right[1], left[2] - right[2]);
}

double norm(const Point3 &vector)
{
    return std::sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
}

} // namespace

ProjectorTable build_projector_table(const PAWSetup &setup, const Structure::LatticeVectors &lattice,
                                     const Atom::AtomicPosition &atom_center, const PeriodicGridPointLocator &locator,
                                     const std::vector<PeriodicNeighbor> &neighbors)
{
    ProjectorTable projector_table;
    projector_table.values =
        ProjectorValueTable(typename ProjectorValueTable::ShapeType{projector_extent(projector_channel_count(setup)),
                                                projector_extent(neighbors.size())});
    projector_table.state_indices.reserve(projector_channel_count(setup));
    projector_table.magnetic_quantum_numbers.reserve(projector_channel_count(setup));
    const Point3 atom_center_cartesian = to_cartesian_point(atom_center);

    for (std::size_t state_index = 0; state_index < setup.projectors_by_state().size(); ++state_index)
    {
        PAWBasisEvaluator evaluator(atom_center, RadialInterpolator(setup.projectors_by_state()[state_index]));
        const int l = setup.states()[state_index].l;

        for (int m = -l; m <= l; ++m)
        {
            const std::size_t channel_index = projector_table.state_indices.size();
            projector_table.state_indices.push_back(state_index);
            projector_table.magnetic_quantum_numbers.push_back(m);

            for (std::size_t neighbor_index = 0; neighbor_index < neighbors.size(); ++neighbor_index)
            {
                const PeriodicNeighbor &neighbor = neighbors[neighbor_index];
                const PeriodicGridPointLocator::Point3 grid_point = locator.point(neighbor.index);
                const Point3 shifted_point = to_shifted_point(make_point3(grid_point[0], grid_point[1], grid_point[2]),
                                                              lattice,
                                                              neighbor.periodic_image);
                const Point3 displacement = subtract_vectors(shifted_point, atom_center_cartesian);
                if (norm(displacement) > setup.cutoff_radius())
                {
                    continue;
                }

                projector_table.values[projector_extent(channel_index), projector_extent(neighbor_index)] =
                    evaluator.evaluate_from_displacement(l, m, displacement);
            }
        }
    }

    return projector_table;
}

DenseMatrix build_projector_basis_overlap_matrix(const DFTGLLHexSpace &space, const PAWSetup &setup, const Atom &atom,
                                                 std::size_t position_index)
{
    mfem::FiniteElementSpace *finite_element_space = const_cast<mfem::FiniteElementSpace *>(&space.space());
    DenseMatrix overlap(
        typename DenseMatrix::ShapeType{projector_extent(projector_channel_count(setup)), projector_extent(space.getDOF())});
    int channel_index = 0;

    for (std::size_t state_index = 0; state_index < setup.projectors_by_state().size(); ++state_index)
    {
        const PAWBasisEvaluator evaluator(atom, position_index, RadialInterpolator(setup.projectors_by_state()[state_index]));
        const int l = setup.states()[state_index].l;

        for (int m = -l; m <= l; ++m)
        {
            mfem::FunctionCoefficient projector_coefficient(
                [evaluator, l, m](const mfem::Vector &point) { return evaluator.evaluate(l, m, to_basis_point(point)); });

            mfem::LinearForm projector_form(finite_element_space);
            projector_form.AddDomainIntegrator(new mfem::DomainLFIntegrator(projector_coefficient));
            projector_form.Assemble();

            mfem::GridFunction projector_grid_function(finite_element_space);
            projector_grid_function = projector_form;

            mfem::Vector true_overlap;
            projector_grid_function.GetTrueDofs(true_overlap);
            for (int dof_index = 0; dof_index < true_overlap.Size(); ++dof_index)
            {
                overlap[channel_index, dof_index] = true_overlap(dof_index);
            }

            ++channel_index;
        }
    }

    return overlap;
}

} // namespace dft
