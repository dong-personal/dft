#include "paw/paw_nonlocal_fixed_operator.hpp"

#include <cmath>

namespace
{

int magnetic_quantum_number_for_state(const dft::PAWState &state)
{
    (void)state;
    return 0;
}

dft::PeriodicGridPointLocator::Vec3 to_locator_point(const dft::Atom::AtomicPosition &point)
{
    return {point[0], point[1], point[2]};
}

} // namespace

PAWNonlocalFixedOperator::PAWNonlocalFixedOperator(const DFTGLLHexSpace &space, const dft::PAWSetup &setup,
                                                   const dft::Atom &atom, std::size_t position_index,
                                                   dft::PeriodicKDTree3D::ImageDepth periodic_images)
    : m_space(space),
      m_setup(setup),
      m_atom_center(to_locator_point(atom.position(position_index))),
      m_locator(space, periodic_images)
{
    m_neighbors = m_locator.RadiusSearch(m_atom_center, m_setup.cutoff_radius());
    m_projector_interpolators.reserve(m_setup.projectors_by_state().size());

    for (const dft::RadialFunction &projector : m_setup.projectors_by_state())
    {
        m_projector_interpolators.emplace_back(projector);
    }

    build_projector_table();
}

void PAWNonlocalFixedOperator::apply(const mfem::Vector &x_true, mfem::Vector &y_true) const
{
    MFEM_VERIFY(x_true.Size() == m_space.getDOF(), "Input size must match true-dof count");
    MFEM_VERIFY(fixed_matrix().size() == m_projector_values.size(),
                "Fixed nonlocal matrix size must match projector count");

    y_true.SetSize(m_space.getDOF());
    y_true = 0.0;

    std::vector<double> projector_coefficients(m_projector_values.size(), 0.0);
    const mfem::Vector &mass_diagonal = m_space.MassDiagTrue();

    for (std::size_t projector_index = 0; projector_index < m_projector_values.size(); ++projector_index)
    {
        double beta_j = 0.0;
        for (std::size_t neighbor_index = 0; neighbor_index < m_neighbors.size(); ++neighbor_index)
        {
            const int point_index = static_cast<int>(m_neighbors[neighbor_index].index);
            beta_j += mass_diagonal(point_index) * m_projector_values[projector_index][neighbor_index] *
                      x_true(point_index);
        }

        projector_coefficients[projector_index] = beta_j;
    }

    std::vector<double> weighted_coefficients(m_projector_values.size(), 0.0);
    for (std::size_t row_index = 0; row_index < m_projector_values.size(); ++row_index)
    {
        double value = 0.0;
        for (std::size_t column_index = 0; column_index < m_projector_values.size(); ++column_index)
        {
            value += fixed_matrix()[row_index][column_index] * projector_coefficients[column_index];
        }

        weighted_coefficients[row_index] = value;
    }

    for (std::size_t neighbor_index = 0; neighbor_index < m_neighbors.size(); ++neighbor_index)
    {
        const int point_index = static_cast<int>(m_neighbors[neighbor_index].index);
        double contribution = 0.0;
        for (std::size_t projector_index = 0; projector_index < m_projector_values.size(); ++projector_index)
        {
            contribution += m_projector_values[projector_index][neighbor_index] *
                            weighted_coefficients[projector_index];
        }

        y_true(point_index) += contribution;
    }
}

PAWNonlocalFixedOperator::Vec3
PAWNonlocalFixedOperator::to_shifted_point(const Vec3 &point, const dft::Structure::LatticeVectors &lattice,
                                           const std::array<int, 3> &periodic_image)
{
    const Vec3 shift = fractional_to_cartesian_shift(lattice, periodic_image);
    return {point[0] + shift[0], point[1] + shift[1], point[2] + shift[2]};
}

PAWNonlocalFixedOperator::Vec3 PAWNonlocalFixedOperator::subtract_vectors(const Vec3 &left, const Vec3 &right)
{
    return {left[0] - right[0], left[1] - right[1], left[2] - right[2]};
}

PAWNonlocalFixedOperator::Vec3
PAWNonlocalFixedOperator::fractional_to_cartesian_shift(const dft::Structure::LatticeVectors &lattice,
                                                        const std::array<int, 3> &periodic_image)
{
    return {
        periodic_image[0] * lattice[0, 0] + periodic_image[1] * lattice[1, 0] + periodic_image[2] * lattice[2, 0],
        periodic_image[0] * lattice[0, 1] + periodic_image[1] * lattice[1, 1] + periodic_image[2] * lattice[2, 1],
        periodic_image[0] * lattice[0, 2] + periodic_image[1] * lattice[1, 2] + periodic_image[2] * lattice[2, 2],
    };
}

double PAWNonlocalFixedOperator::norm(const Vec3 &vector)
{
    return std::sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
}

void PAWNonlocalFixedOperator::build_projector_table()
{
    m_projector_values.assign(m_projector_interpolators.size(), std::vector<double>(m_neighbors.size(), 0.0));

    const auto &points = m_locator.points();
    const auto &lattice = m_space.MeshSource().lattice();

    for (std::size_t state_index = 0; state_index < m_projector_interpolators.size(); ++state_index)
    {
        dft::Atom::AtomicPosition evaluator_center(typename dft::Atom::AtomicPosition::ShapeType{3});
        evaluator_center[0] = m_atom_center[0];
        evaluator_center[1] = m_atom_center[1];
        evaluator_center[2] = m_atom_center[2];

        dft::PAWBasisEvaluator evaluator(evaluator_center, m_projector_interpolators[state_index]);
        const int l = m_setup.states()[state_index].l;
        const int m = magnetic_quantum_number_for_state(m_setup.states()[state_index]);

        for (std::size_t neighbor_index = 0; neighbor_index < m_neighbors.size(); ++neighbor_index)
        {
            const dft::PeriodicNeighbor &neighbor = m_neighbors[neighbor_index];
            const Vec3 shifted_point = to_shifted_point(points[neighbor.index], lattice, neighbor.periodic_image);
            const Vec3 displacement = subtract_vectors(shifted_point, m_atom_center);
            if (norm(displacement) > m_setup.cutoff_radius())
            {
                continue;
            }

            m_projector_values[state_index][neighbor_index] = evaluator.evaluate_from_displacement(l, m, displacement);
        }
    }
}
