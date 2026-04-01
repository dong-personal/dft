#include "paw/interpolator.hpp"

#include "paw/spherical_harmonics.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace dft
{

namespace
{

constexpr double kZeroRadiusTolerance = 0.0;

PAWBasisEvaluator::Point3 make_point3(double x, double y, double z)
{
    PAWBasisEvaluator::Point3 point(typename PAWBasisEvaluator::Point3::ShapeType{3});
    point[0] = x;
    point[1] = y;
    point[2] = z;
    return point;
}

PAWBasisEvaluator::Point3 subtract_vectors(const PAWBasisEvaluator::Point3 &left,
                                           const PAWBasisEvaluator::Point3 &right)
{
    return make_point3(left[0] - right[0], left[1] - right[1], left[2] - right[2]);
}

double norm(const PAWBasisEvaluator::Point3 &vector)
{
    return std::sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
}

} // namespace

RadialInterpolator::RadialInterpolator(const RadialFunction &radial_function)
    : m_radii(radial_function.radii), m_values(radial_function.values)
{
    validate();
}

RadialInterpolator::RadialInterpolator(std::vector<double> radii, std::vector<double> values)
    : m_radii(std::move(radii)), m_values(std::move(values))
{
    validate();
}

double RadialInterpolator::evaluate(double radius) const
{
    if (radius <= m_radii.front())
    {
        return m_values.front();
    }

    if (radius >= m_radii.back())
    {
        return m_values.back();
    }

    const auto upper_iterator = std::upper_bound(m_radii.begin(), m_radii.end(), radius);
    const std::size_t upper_index = static_cast<std::size_t>(upper_iterator - m_radii.begin());
    const std::size_t lower_index = upper_index - 1;

    const double r0 = m_radii[lower_index];
    const double r1 = m_radii[upper_index];
    const double v0 = m_values[lower_index];
    const double v1 = m_values[upper_index];
    const double interpolation_fraction = (radius - r0) / (r1 - r0);

    return (1.0 - interpolation_fraction) * v0 + interpolation_fraction * v1;
}

void RadialInterpolator::validate() const
{
    if (m_radii.empty() || m_values.empty())
    {
        throw std::runtime_error("RadialInterpolator requires non-empty radial data");
    }

    if (m_radii.size() != m_values.size())
    {
        throw std::runtime_error("RadialInterpolator radial grid/value size mismatch");
    }

    if (m_radii.size() < 2)
    {
        throw std::runtime_error("RadialInterpolator requires at least two grid points");
    }

    for (std::size_t index = 1; index < m_radii.size(); ++index)
    {
        if (!(m_radii[index] > m_radii[index - 1]))
        {
            throw std::runtime_error("RadialInterpolator requires a strictly increasing radial grid");
        }
    }
}

PAWBasisEvaluator::PAWBasisEvaluator(Atom::AtomicPosition center, RadialInterpolator radial_interpolator)
    : m_center(std::move(center)), m_radial_interpolator(std::move(radial_interpolator))
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

double PAWBasisEvaluator::evaluate(int l, int m, const Point3 &point) const
{
    return evaluate_from_displacement(l, m, subtract_vectors(point, m_center));
}

double PAWBasisEvaluator::evaluate_from_displacement(int l, int m, const Point3 &displacement) const
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
