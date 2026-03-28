#pragma once

#include "atom.h"
#include "radial_interpolator.hpp"

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

    double Evaluate(int l, int m, const Vec3 &point) const;
    double EvaluateFromDisplacement(int l, int m, const Vec3 &displacement) const;

    const Vec3 &center() const { return center_; }
    const RadialInterpolator &radial_interpolator() const { return radial_interpolator_; }

  private:
    Vec3 center_;
    RadialInterpolator radial_interpolator_;
};

} // namespace dft
