#include "paw/spherical_harmonics.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dft
{

namespace
{

constexpr double kPi = 3.141592653589793238462643383279502884;

double double_factorial(int n)
{
    if (n <= 0)
    {
        return 1.0;
    }

    double result = 1.0;
    for (int value = n; value > 0; value -= 2)
    {
        result *= static_cast<double>(value);
    }

    return result;
}

double normalization_constant(int l, int absolute_m)
{
    const double numerator = (2.0 * static_cast<double>(l) + 1.0) * factorial(l - absolute_m);
    const double denominator = 4.0 * kPi * factorial(l + absolute_m);
    return std::sqrt(numerator / denominator);
}

} // namespace

double factorial(int n)
{
    if (n < 0)
    {
        throw std::invalid_argument("Factorial is undefined for negative integers");
    }

    double result = 1.0;
    for (int value = 2; value <= n; ++value)
    {
        result *= static_cast<double>(value);
    }

    return result;
}

std::pair<double, double> cartesian_to_spherical_angles(double x, double y, double z)
{
    const double radius = std::sqrt(x * x + y * y + z * z);
    if (radius == 0.0)
    {
        throw std::invalid_argument("Spherical angles are undefined at the origin");
    }

    const double cos_theta = std::clamp(z / radius, -1.0, 1.0);
    const double theta = std::acos(cos_theta);
    const double phi = std::atan2(y, x);

    return {theta, phi};
}

double associated_legendre_polynomial(int l, int m, double x)
{
    if (l < 0)
    {
        throw std::invalid_argument("Associated Legendre polynomial requires l >= 0");
    }

    if (std::abs(m) > l)
    {
        throw std::invalid_argument("Associated Legendre polynomial requires |m| <= l");
    }

    if (x < -1.0 || x > 1.0)
    {
        throw std::invalid_argument("Associated Legendre polynomial requires x in [-1, 1]");
    }

    const int absolute_m = std::abs(m);

    double p_mm = 1.0;
    if (absolute_m > 0)
    {
        const double sin_theta = std::sqrt(std::max(0.0, 1.0 - x * x));
        p_mm = std::pow(-1.0, absolute_m) * double_factorial(2 * absolute_m - 1) * std::pow(sin_theta, absolute_m);
    }

    if (l == absolute_m)
    {
        if (m >= 0)
        {
            return p_mm;
        }

        const double sign = (absolute_m % 2 == 0) ? 1.0 : -1.0;
        return sign * factorial(l - absolute_m) / factorial(l + absolute_m) * p_mm;
    }

    const double p_mmp1 = x * (2 * absolute_m + 1) * p_mm;
    if (l == absolute_m + 1)
    {
        if (m >= 0)
        {
            return p_mmp1;
        }

        const double sign = (absolute_m % 2 == 0) ? 1.0 : -1.0;
        return sign * factorial(l - absolute_m) / factorial(l + absolute_m) * p_mmp1;
    }

    double p_lm_minus_two = p_mm;
    double p_lm_minus_one = p_mmp1;
    double p_lm = 0.0;
    for (int ell = absolute_m + 2; ell <= l; ++ell)
    {
        p_lm = ((2.0 * ell - 1.0) * x * p_lm_minus_one - (ell + absolute_m - 1.0) * p_lm_minus_two) /
               (ell - absolute_m);
        p_lm_minus_two = p_lm_minus_one;
        p_lm_minus_one = p_lm;
    }

    if (m >= 0)
    {
        return p_lm;
    }

    const double sign = (absolute_m % 2 == 0) ? 1.0 : -1.0;
    return sign * factorial(l - absolute_m) / factorial(l + absolute_m) * p_lm;
}

double evaluate_real_spherical_harmonic_from_angles(int l, int m, double theta, double phi)
{
    if (l < 0)
    {
        throw std::invalid_argument("Real spherical harmonics require l >= 0");
    }

    if (std::abs(m) > l)
    {
        throw std::invalid_argument("Real spherical harmonics require |m| <= l");
    }

    const int absolute_m = std::abs(m);
    const double legendre = associated_legendre_polynomial(l, absolute_m, std::cos(theta));
    const double normalization = normalization_constant(l, absolute_m);

    if (m == 0)
    {
        return normalization * legendre;
    }

    if (m > 0)
    {
        return std::sqrt(2.0) * normalization * legendre * std::cos(static_cast<double>(absolute_m) * phi);
    }

    return std::sqrt(2.0) * normalization * legendre * std::sin(static_cast<double>(absolute_m) * phi);
}

double evaluate_real_spherical_harmonic(int l, int m, double x, double y, double z)
{
    const auto [theta, phi] = cartesian_to_spherical_angles(x, y, z);
    return evaluate_real_spherical_harmonic_from_angles(l, m, theta, phi);
}

} // namespace dft
