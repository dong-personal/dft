#include "sphericalHarmonics.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dft
{

namespace
{

constexpr double PI = 3.141592653589793238462643383279502884;

double DoubleFactorial(int n)
{
    if (n <= 0)
    {
        return 1.0;
    }

    double result = 1.0;
    for (int k = n; k > 0; k -= 2)
    {
        result *= static_cast<double>(k);
    }
    return result;
}

double NormalizationConstant(int l, int m_abs)
{
    const double numerator = (2.0 * static_cast<double>(l) + 1.0) * Factorial(l - m_abs);
    const double denominator = 4.0 * PI * Factorial(l + m_abs);
    return std::sqrt(numerator / denominator);
}

} // namespace

double Factorial(int n)
{
    if (n < 0)
    {
        throw std::invalid_argument("Factorial is undefined for negative integers");
    }

    double result = 1.0;
    for (int k = 2; k <= n; ++k)
    {
        result *= static_cast<double>(k);
    }
    return result;
}

std::pair<double, double> CartesianToSphericalAngles(double x, double y, double z)
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

double AssociatedLegendrePolynomial(int l, int m, double x)
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

    const int m_abs = std::abs(m);

    double p_mm = 1.0;
    if (m_abs > 0)
    {
        const double sin_theta = std::sqrt(std::max(0.0, 1.0 - x * x));
        p_mm = std::pow(-1.0, m_abs) * DoubleFactorial(2 * m_abs - 1) * std::pow(sin_theta, m_abs);
    }

    if (l == m_abs)
    {
        if (m >= 0)
        {
            return p_mm;
        }
        const double sign = (m_abs % 2 == 0) ? 1.0 : -1.0;
        return sign * Factorial(l - m_abs) / Factorial(l + m_abs) * p_mm;
    }

    double p_mmp1 = x * (2 * m_abs + 1) * p_mm;
    if (l == m_abs + 1)
    {
        if (m >= 0)
        {
            return p_mmp1;
        }
        const double sign = (m_abs % 2 == 0) ? 1.0 : -1.0;
        return sign * Factorial(l - m_abs) / Factorial(l + m_abs) * p_mmp1;
    }

    double p_lm_minus2 = p_mm;
    double p_lm_minus1 = p_mmp1;
    double p_lm = 0.0;
    for (int ell = m_abs + 2; ell <= l; ++ell)
    {
        p_lm = ((2.0 * ell - 1.0) * x * p_lm_minus1 - (ell + m_abs - 1.0) * p_lm_minus2) /
               (ell - m_abs);
        p_lm_minus2 = p_lm_minus1;
        p_lm_minus1 = p_lm;
    }

    if (m >= 0)
    {
        return p_lm;
    }

    const double sign = (m_abs % 2 == 0) ? 1.0 : -1.0;
    return sign * Factorial(l - m_abs) / Factorial(l + m_abs) * p_lm;
}

double EvaluateRealSphericalHarmonicFromAngles(int l, int m, double theta, double phi)
{
    if (l < 0)
    {
        throw std::invalid_argument("Real spherical harmonics require l >= 0");
    }
    if (std::abs(m) > l)
    {
        throw std::invalid_argument("Real spherical harmonics require |m| <= l");
    }

    const int m_abs = std::abs(m);
    const double legendre = AssociatedLegendrePolynomial(l, m_abs, std::cos(theta));
    const double normalization = NormalizationConstant(l, m_abs);

    if (m == 0)
    {
        return normalization * legendre;
    }
    if (m > 0)
    {
        return std::sqrt(2.0) * normalization * legendre * std::cos(static_cast<double>(m_abs) * phi);
    }
    return std::sqrt(2.0) * normalization * legendre * std::sin(static_cast<double>(m_abs) * phi);
}

double EvaluateRealSphericalHarmonic(int l, int m, double x, double y, double z)
{
    const auto [theta, phi] = CartesianToSphericalAngles(x, y, z);
    return EvaluateRealSphericalHarmonicFromAngles(l, m, theta, phi);
}

} // namespace dft
