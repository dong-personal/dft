#pragma once
#ifndef DFT_INTERPOLATOR_HPP
#define DFT_INTERPOLATOR_HPP
#include "atom.h"
#include "paw/paw_setup.hpp"

#include <cstddef>
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

class PAWBasisEvaluator
{
  public:
    using Point3 = Atom::AtomicPosition;

    PAWBasisEvaluator(Atom::AtomicPosition center, RadialInterpolator radial_interpolator);
    PAWBasisEvaluator(const Atom &atom, std::size_t position_index, RadialInterpolator radial_interpolator);
    PAWBasisEvaluator(const Atom &atom, RadialInterpolator radial_interpolator);

    double evaluate(int l, int m, const Point3 &point) const;
    double evaluate_from_displacement(int l, int m, const Point3 &displacement) const;

    const Point3 &center() const { return m_center; }
    const RadialInterpolator &radial_interpolator() const { return m_radial_interpolator; }

  private:
    Point3 m_center{typename Point3::ShapeType{3}};
    RadialInterpolator m_radial_interpolator;
};

} // namespace dft

#endif // DFT_INTERPOLATOR_HPP
