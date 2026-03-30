#pragma once

#include "atom.h"
#include "paw/radial_interpolator.hpp"

#include <array>
#include <cstddef>

namespace dft
{

class PAWBasisEvaluator
{
  public:
    using Vec3 = std::array<double, 3>;

    PAWBasisEvaluator(Atom::AtomicPosition center, RadialInterpolator radial_interpolator);
    PAWBasisEvaluator(const Atom &atom, std::size_t position_index, RadialInterpolator radial_interpolator);
    PAWBasisEvaluator(const Atom &atom, RadialInterpolator radial_interpolator);

    double evaluate(int l, int m, const Vec3 &point) const;
    double evaluate_from_displacement(int l, int m, const Vec3 &displacement) const;

    const Vec3 &center() const { return m_center; }
    const RadialInterpolator &radial_interpolator() const { return m_radial_interpolator; }

  private:
    Vec3 m_center{};
    RadialInterpolator m_radial_interpolator;
};

} // namespace dft
