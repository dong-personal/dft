#include "paw_basis_evaluator.hpp"

#include "spherical_harmonics.hpp"

#include <cmath>
#include <utility>

namespace dft
{

namespace
{

PAWBasisEvaluator::Vec3 Subtract(const PAWBasisEvaluator::Vec3 &a, const PAWBasisEvaluator::Vec3 &b)
{
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

double Norm(const PAWBasisEvaluator::Vec3 &v)
{
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

} // namespace

PAWBasisEvaluator::PAWBasisEvaluator(Atom::AtomicPosition center, RadialInterpolator radial_interpolator)
    : center_{center[0], center[1], center[2]}, radial_interpolator_(std::move(radial_interpolator))
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

double PAWBasisEvaluator::Evaluate(int l, int m, const Vec3 &point) const
{
    return EvaluateFromDisplacement(l, m, Subtract(point, center_));
}

double PAWBasisEvaluator::EvaluateFromDisplacement(int l, int m, const Vec3 &displacement) const
{
    const double radius = Norm(displacement);
    if (radius == 0.0)
    {
        if (l == 0 && m == 0)
        {
            return radial_interpolator_.Evaluate(0.0) * EvaluateRealSphericalHarmonicFromAngles(0, 0, 0.0, 0.0);
        }
        return 0.0;
    }

    if (radius > radial_interpolator_.r().back())
    {
        return 0.0;
    }

    const double radial_value = radial_interpolator_.Evaluate(radius);
    const double angular_value =
        EvaluateRealSphericalHarmonic(l, m, displacement[0], displacement[1], displacement[2]);

    return radial_value * angular_value;
}

} // namespace dft
