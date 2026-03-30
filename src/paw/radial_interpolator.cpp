#include "paw/radial_interpolator.hpp"

#include <algorithm>
#include <stdexcept>

namespace dft
{

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

} // namespace dft
