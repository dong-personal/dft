#pragma once
#ifndef DFT_SPHERICAL_HARMONICS_HPP
#define DFT_SPHERICAL_HARMONICS_HPP
#include <utility>

namespace dft
{

double factorial(int n);
std::pair<double, double> cartesian_to_spherical_angles(double x, double y, double z);

double associated_legendre_polynomial(int l, int m, double x);
double evaluate_real_spherical_harmonic_from_angles(int l, int m, double theta, double phi);
double evaluate_real_spherical_harmonic(int l, int m, double x, double y, double z);

} // namespace dft
#endif // DFT_SPHERICAL_HARMONICS_HPP
