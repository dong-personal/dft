#include "paw/paw_basis_evaluator.hpp"

#include "paw/spherical_harmonics.hpp"

#include <cmath>
#include <utility>

namespace dft
{

namespace
{

constexpr double kZeroRadiusTolerance = 0.0;

PAWBasisEvaluator::Vec3 subtract_vectors(const PAWBasisEvaluator::Vec3 &left, const PAWBasisEvaluator::Vec3 &right)
{
    return {left[0] - right[0], left[1] - right[1], left[2] - right[2]};
}

double norm(const PAWBasisEvaluator::Vec3 &vector)
{
    return std::sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
}

} // namespace

PAWBasisEvaluator::PAWBasisEvaluator(Atom::AtomicPosition center, RadialInterpolator radial_interpolator)
    : m_center{center[0], center[1], center[2]}, m_radial_interpolator(std::move(radial_interpolator))
{
}

PAWBasisEvaluator::PAWBasisEvaluator(const Atom &atom, std::size_t position_index,
                                     RadialInterpolator radial_interpolator)
    : PAWBasisEvaluator(atom.position(position_index), std::move(radial_interpolator))
{
}

PAWBasisEvaluator::PAWBasisEvaluator(const Atom &atom, RadialInterpolator radial_interpolator)
    : PAWBasisEvaluator(atom, 0, std::move(radial_interpolator))
{
}

double PAWBasisEvaluator::evaluate(int l, int m, const Vec3 &point) const
{
    return evaluate_from_displacement(l, m, subtract_vectors(point, m_center));
}

double PAWBasisEvaluator::evaluate_from_displacement(int l, int m, const Vec3 &displacement) const
{
    const double radius = norm(displacement);
    if (radius == kZeroRadiusTolerance)
    {
        if (l == 0 && m == 0)
        {
            return m_radial_interpolator.evaluate(0.0) *
                   evaluate_real_spherical_harmonic_from_angles(0, 0, 0.0, 0.0);
        }

        return 0.0;
    }

    if (radius > m_radial_interpolator.radii().back())
    {
        return 0.0;
    }

    const double radial_value = m_radial_interpolator.evaluate(radius);
    const double angular_value =
        evaluate_real_spherical_harmonic(l, m, displacement[0], displacement[1], displacement[2]);

    return radial_value * angular_value;
}

} // namespace dft
