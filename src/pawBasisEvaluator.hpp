#pragma once

#include "Atom.h"
#include "radialInterpolator.hpp"

#include <array>

namespace dft
{

class PAWBasisEvaluator
{
  public:
    using Vec3 = std::array<double, 3>;

    PAWBasisEvaluator(Atom::Vec3 center, RadialInterpolator radial_interpolator);
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
