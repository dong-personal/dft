#pragma once

#include <utility>

namespace dft
{

double Factorial(int n);
std::pair<double, double> CartesianToSphericalAngles(double x, double y, double z);

double AssociatedLegendrePolynomial(int l, int m, double x);
double EvaluateRealSphericalHarmonicFromAngles(int l, int m, double theta, double phi);
double EvaluateRealSphericalHarmonic(int l, int m, double x, double y, double z);

} // namespace dft
