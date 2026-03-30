#pragma once

#include "paw/paw.hpp"

#include <vector>

namespace dft
{

class RadialInterpolator
{
  public:
    explicit RadialInterpolator(const RadialFunction &radial_function);
    RadialInterpolator(std::vector<double> radii, std::vector<double> values);

    double evaluate(double radius) const;

    const std::vector<double> &radii() const { return m_radii; }
    const std::vector<double> &values() const { return m_values; }

  private:
    void validate() const;

    std::vector<double> m_radii;
    std::vector<double> m_values;
};

} // namespace dft
