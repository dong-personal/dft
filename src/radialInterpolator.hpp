#pragma once

#include "PAW.h"

#include <vector>

namespace dft
{

class RadialInterpolator
{
  public:
    explicit RadialInterpolator(const RadialFunction &radial_function);
    RadialInterpolator(std::vector<double> r, std::vector<double> values);

    double Evaluate(double radius) const;

    const std::vector<double> &r() const { return r_; }
    const std::vector<double> &values() const { return values_; }

  private:
    void Validate_() const;

    std::vector<double> r_;
    std::vector<double> values_;
};

} // namespace dft
