#include "sphericalHarmonics.hpp"

#include <cmath>
#include <iostream>

namespace
{

bool NearlyEqual(double a, double b, double tol = 1e-12)
{
    return std::abs(a - b) <= tol;
}

constexpr double PI = 3.141592653589793238462643383279502884;

} // namespace

int main()
{
    const double y00 = dft::EvaluateRealSphericalHarmonicFromAngles(0, 0, PI / 2.0, 0.0);
    const double expected_y00 = 0.5 / std::sqrt(PI);

    const double y10_z = dft::EvaluateRealSphericalHarmonic(1, 0, 0.0, 0.0, 1.0);
    const double expected_y10_z = std::sqrt(3.0 / (4.0 * PI));

    const double y11_x = dft::EvaluateRealSphericalHarmonic(1, 1, 1.0, 0.0, 0.0);
    const double expected_y11_x = -std::sqrt(3.0 / (4.0 * PI));

    const double y1m1_y = dft::EvaluateRealSphericalHarmonic(1, -1, 0.0, 1.0, 0.0);
    const double expected_y1m1_y = -std::sqrt(3.0 / (4.0 * PI));

    const double y20_z = dft::EvaluateRealSphericalHarmonic(2, 0, 0.0, 0.0, 1.0);
    const double expected_y20_z = std::sqrt(5.0 / (4.0 * PI));

    std::cout << "Y00 = " << y00 << '\n';
    std::cout << "Y10(z-axis) = " << y10_z << '\n';
    std::cout << "Y11(x-axis) = " << y11_x << '\n';
    std::cout << "Y1-1(y-axis) = " << y1m1_y << '\n';
    std::cout << "Y20(z-axis) = " << y20_z << '\n';

    if (!NearlyEqual(y00, expected_y00))
    {
        std::cerr << "Y00 normalization check failed" << std::endl;
        return 1;
    }
    if (!NearlyEqual(y10_z, expected_y10_z))
    {
        std::cerr << "Y10 z-axis check failed" << std::endl;
        return 1;
    }
    if (!NearlyEqual(y11_x, expected_y11_x))
    {
        std::cerr << "Y11 x-axis check failed" << std::endl;
        return 1;
    }
    if (!NearlyEqual(y1m1_y, expected_y1m1_y))
    {
        std::cerr << "Y1-1 y-axis check failed" << std::endl;
        return 1;
    }
    if (!NearlyEqual(y20_z, expected_y20_z))
    {
        std::cerr << "Y20 z-axis check failed" << std::endl;
        return 1;
    }

    return 0;
}
