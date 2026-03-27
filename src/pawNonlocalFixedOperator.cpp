#include "pawNonlocalFixedOperator.hpp"

#include "pawBasisEvaluator.hpp"

#include <cmath>
#include <utility>

namespace
{

int MagneticQuantumNumberForState(const dft::PAWState &state)
{
    (void)state;
    return 0;
}

} // namespace

PAWNonlocalFixedOperator::PAWNonlocalFixedOperator(const DFTGLLHexSpace &space, const dft::PAWSetup &setup,
                                                   const dft::Atom &atom,
                                                   dft::PeriodicKDTree3D::ImageDepth periodic_images)
    : space_(space), setup_(setup), atom_(atom), locator_(space, periodic_images)
{
    neighbors_ = locator_.RadiusSearch(atom.position(), setup_.cutoff_radius());
    projector_interpolators_.reserve(setup_.projectors_by_state().size());
    for (const auto &projector : setup_.projectors_by_state())
    {
        projector_interpolators_.emplace_back(projector);
    }
    BuildProjectorTable_();
}

void PAWNonlocalFixedOperator::Apply(const mfem::Vector &x_true, mfem::Vector &y_true) const
{
    MFEM_VERIFY(x_true.Size() == space_.getDOF(), "Input size must match true-dof count");
    MFEM_VERIFY(fixed_matrix().size() == projector_values_.size(),
                "Fixed nonlocal matrix size must match projector count");

    y_true.SetSize(space_.getDOF());
    y_true = 0.0;

    std::vector<double> projector_coeffs(projector_values_.size(), 0.0);
    const mfem::Vector &mass_diag = space_.MassDiagTrue();

    for (std::size_t j = 0; j < projector_values_.size(); ++j)
    {
        double beta_j = 0.0;
        for (std::size_t hit = 0; hit < neighbors_.size(); ++hit)
        {
            const int idx = static_cast<int>(neighbors_[hit].index);
            beta_j += mass_diag(idx) * projector_values_[j][hit] * x_true(idx);
        }
        projector_coeffs[j] = beta_j;
    }

    std::vector<double> weighted_coeffs(projector_values_.size(), 0.0);
    for (std::size_t i = 0; i < projector_values_.size(); ++i)
    {
        double value = 0.0;
        for (std::size_t j = 0; j < projector_values_.size(); ++j)
        {
            value += fixed_matrix()[i][j] * projector_coeffs[j];
        }
        weighted_coeffs[i] = value;
    }

    for (std::size_t hit = 0; hit < neighbors_.size(); ++hit)
    {
        const int idx = static_cast<int>(neighbors_[hit].index);
        double contribution = 0.0;
        for (std::size_t i = 0; i < projector_values_.size(); ++i)
        {
            contribution += projector_values_[i][hit] * weighted_coeffs[i];
        }
        y_true(idx) += contribution;
    }
}

PAWNonlocalFixedOperator::Vec3 PAWNonlocalFixedOperator::ToShiftedPoint_(const Vec3 &point,
                                                                         const dft::Structure::Mat3 &lattice,
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
PAWNonlocalFixedOperator::FractionalToCartesianShift_(const dft::Structure::Mat3 &lattice,
                                                      const std::array<int, 3> &periodic_image)
{
    return {
        periodic_image[0] * lattice[0][0] + periodic_image[1] * lattice[1][0] + periodic_image[2] * lattice[2][0],
        periodic_image[0] * lattice[0][1] + periodic_image[1] * lattice[1][1] + periodic_image[2] * lattice[2][1],
        periodic_image[0] * lattice[0][2] + periodic_image[1] * lattice[1][2] + periodic_image[2] * lattice[2][2],
    };
}

double PAWNonlocalFixedOperator::Norm_(const Vec3 &v)
{
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

void PAWNonlocalFixedOperator::BuildProjectorTable_()
{
    projector_values_.assign(projector_interpolators_.size(), std::vector<double>(neighbors_.size(), 0.0));

    const auto &points = locator_.points();
    const auto &lattice = space_.MeshSource().lattice();
    const Vec3 atom_center = atom_.position();

    for (std::size_t state_index = 0; state_index < projector_interpolators_.size(); ++state_index)
    {
        dft::PAWBasisEvaluator evaluator(atom_center, projector_interpolators_[state_index]);
        const int l = setup_.states()[state_index].l;
        const int m = MagneticQuantumNumberForState(setup_.states()[state_index]);

        for (std::size_t hit = 0; hit < neighbors_.size(); ++hit)
        {
            const auto &neighbor = neighbors_[hit];
            const Vec3 shifted_point =
                ToShiftedPoint_(points[neighbor.index], lattice, neighbor.periodic_image);
            const Vec3 displacement = Subtract_(shifted_point, atom_center);
            if (Norm_(displacement) > setup_.cutoff_radius())
            {
                continue;
            }
            projector_values_[state_index][hit] = evaluator.EvaluateFromDisplacement(l, m, displacement);
        }
    }
}
