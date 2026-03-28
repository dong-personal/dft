#include "paw_nonlocal_fixed_operator.hpp"

#include "paw_basis_evaluator.hpp"

#include <cmath>

namespace
{

int MagneticQuantumNumberForState(const dft::PAWState &state)
{
    (void)state;
    return 0;
}

dft::PeriodicGridPointLocator::Vec3 ToLocatorPoint(const dft::Atom::AtomicPosition &point)
{
    return {point[0], point[1], point[2]};
}

} // namespace

PAWNonlocalFixedOperator::PAWNonlocalFixedOperator(const DFTGLLHexSpace &space, const dft::PAWSetup &setup,
                                                   const dft::Atom &atom, std::size_t position_index,
                                                   dft::PeriodicKDTree3D::ImageDepth periodic_images)
    : m_space(space),
      m_setup(setup),
      m_atom_center(ToLocatorPoint(atom.position(position_index))),
      m_locator(space, periodic_images)
{
    m_neighbors = m_locator.RadiusSearch(m_atom_center, m_setup.cutoff_radius());
    m_projector_interpolators.reserve(m_setup.projectors_by_state().size());
    for (const auto &projector : m_setup.projectors_by_state())
    {
        m_projector_interpolators.emplace_back(projector);
    }
    BuildProjectorTable_();
}

void PAWNonlocalFixedOperator::Apply(const mfem::Vector &x_true, mfem::Vector &y_true) const
{
    MFEM_VERIFY(x_true.Size() == m_space.getDOF(), "Input size must match true-dof count");
    MFEM_VERIFY(fixed_matrix().size() == m_projector_values.size(),
                "Fixed nonlocal matrix size must match projector count");

    y_true.SetSize(m_space.getDOF());
    y_true = 0.0;

    std::vector<double> projector_coeffs(m_projector_values.size(), 0.0);
    const mfem::Vector &mass_diag = m_space.MassDiagTrue();

    for (std::size_t j = 0; j < m_projector_values.size(); ++j)
    {
        double beta_j = 0.0;
        for (std::size_t hit = 0; hit < m_neighbors.size(); ++hit)
        {
            const int idx = static_cast<int>(m_neighbors[hit].index);
            beta_j += mass_diag(idx) * m_projector_values[j][hit] * x_true(idx);
        }
        projector_coeffs[j] = beta_j;
    }

    std::vector<double> weighted_coeffs(m_projector_values.size(), 0.0);
    for (std::size_t i = 0; i < m_projector_values.size(); ++i)
    {
        double value = 0.0;
        for (std::size_t j = 0; j < m_projector_values.size(); ++j)
        {
            value += fixed_matrix()[i][j] * projector_coeffs[j];
        }
        weighted_coeffs[i] = value;
    }

    for (std::size_t hit = 0; hit < m_neighbors.size(); ++hit)
    {
        const int idx = static_cast<int>(m_neighbors[hit].index);
        double contribution = 0.0;
        for (std::size_t i = 0; i < m_projector_values.size(); ++i)
        {
            contribution += m_projector_values[i][hit] * weighted_coeffs[i];
        }
        y_true(idx) += contribution;
    }
}

PAWNonlocalFixedOperator::Vec3 PAWNonlocalFixedOperator::ToShiftedPoint_(const Vec3 &point,
                                                                         const dft::Structure::LatticeVectors &lattice,
                                                                         const std::array<int, 3> &periodic_image)
{
    const Vec3 shift = FractionalToCartesianShift_(lattice, periodic_image);
    return {point[0] + shift[0], point[1] + shift[1], point[2] + shift[2]};
}

PAWNonlocalFixedOperator::Vec3 PAWNonlocalFixedOperator::Subtract_(const Vec3 &a, const Vec3 &b)
{
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

PAWNonlocalFixedOperator::Vec3
PAWNonlocalFixedOperator::FractionalToCartesianShift_(const dft::Structure::LatticeVectors &lattice,
                                                      const std::array<int, 3> &periodic_image)
{
    return {
        periodic_image[0] * lattice[0, 0] + periodic_image[1] * lattice[1, 0] + periodic_image[2] * lattice[2, 0],
        periodic_image[0] * lattice[0, 1] + periodic_image[1] * lattice[1, 1] + periodic_image[2] * lattice[2, 1],
        periodic_image[0] * lattice[0, 2] + periodic_image[1] * lattice[1, 2] + periodic_image[2] * lattice[2, 2],
    };
}

double PAWNonlocalFixedOperator::Norm_(const Vec3 &v)
{
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

void PAWNonlocalFixedOperator::BuildProjectorTable_()
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
        const int m = MagneticQuantumNumberForState(m_setup.states()[state_index]);

        for (std::size_t hit = 0; hit < m_neighbors.size(); ++hit)
        {
            const auto &neighbor = m_neighbors[hit];
            const Vec3 shifted_point =
                ToShiftedPoint_(points[neighbor.index], lattice, neighbor.periodic_image);
            const Vec3 displacement = Subtract_(shifted_point, m_atom_center);
            if (Norm_(displacement) > m_setup.cutoff_radius())
            {
                continue;
            }
            m_projector_values[state_index][hit] = evaluator.EvaluateFromDisplacement(l, m, displacement);
        }
    }
}
