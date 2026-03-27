#include "radialInterpolator.hpp"

#include <algorithm>
#include <stdexcept>

namespace dft
{

RadialInterpolator::RadialInterpolator(const RadialFunction &radial_function)
    : r_(radial_function.r), values_(radial_function.values)
{
    Validate_();
}

RadialInterpolator::RadialInterpolator(std::vector<double> r, std::vector<double> values)
    : r_(std::move(r)), values_(std::move(values))
{
    Validate_();
}

double RadialInterpolator::Evaluate(double radius) const
{
    if (radius <= r_.front())
    {
        return values_.front();
    }
    if (radius >= r_.back())
    {
        return values_.back();
    }

    auto upper = std::upper_bound(r_.begin(), r_.end(), radius);
    const std::size_t hi = static_cast<std::size_t>(upper - r_.begin());
    const std::size_t lo = hi - 1;

    const double r0 = r_[lo];
    const double r1 = r_[hi];
    const double v0 = values_[lo];
    const double v1 = values_[hi];

    const double t = (radius - r0) / (r1 - r0);
    return (1.0 - t) * v0 + t * v1;
}

void RadialInterpolator::Validate_() const
{
    if (r_.empty() || values_.empty())
    {
        throw std::runtime_error("RadialInterpolator requires non-empty radial data");
    }
    if (r_.size() != values_.size())
    {
        throw std::runtime_error("RadialInterpolator radial grid/value size mismatch");
    }
    if (r_.size() < 2)
    {
        throw std::runtime_error("RadialInterpolator requires at least two grid points");
    }
    for (std::size_t i = 1; i < r_.size(); ++i)
    {
        if (!(r_[i] > r_[i - 1]))
        {
            throw std::runtime_error("RadialInterpolator requires a strictly increasing radial grid");
        }
    }
}

} // namespace dft
