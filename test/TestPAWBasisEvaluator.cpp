#include "paw/paw_basis_evaluator.hpp"

#include <cmath>
#include <iostream>
#include <vector>

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
    const std::vector<double> radial_grid{0.0, 1.0, 2.0};
    const std::vector<double> radial_values{2.0, 2.0, 2.0};
    dft::RadialInterpolator radial(radial_grid, radial_values);

    dft::Atom::AtomicPosition atom_position(typename dft::Atom::AtomicPosition::ShapeType{3});
    atom_position[0] = 0.0;
    atom_position[1] = 0.0;
    atom_position[2] = 0.0;

    dft::Atom atom("X", 1);
    atom.set_position(0, atom_position);
    const dft::PAWBasisEvaluator evaluator(atom, radial);

    const double expected_y00 = 1.0 / (2.0 * std::sqrt(PI));
    const double expected_y10 = std::sqrt(3.0 / (4.0 * PI));
    const double expected_y11 = -std::sqrt(3.0 / (4.0 * PI));

    const double s_origin = evaluator.evaluate(0, 0, {0.0, 0.0, 0.0});
    const double pz_axis = evaluator.evaluate(1, 0, {0.0, 0.0, 1.0});
    const double px_axis = evaluator.evaluate(1, 1, {1.0, 0.0, 0.0});
    const double outside_cutoff = evaluator.evaluate(0, 0, {0.0, 0.0, 3.0});

    std::cout << "s(origin) = " << s_origin << '\n';
    std::cout << "p_z(z-axis) = " << pz_axis << '\n';
    std::cout << "p_x(x-axis) = " << px_axis << '\n';
    std::cout << "outside cutoff = " << outside_cutoff << '\n';

    if (!NearlyEqual(s_origin, 2.0 * expected_y00))
    {
        std::cerr << "s-like basis value at the origin is incorrect" << std::endl;
        return 1;
    }
    if (!NearlyEqual(pz_axis, 2.0 * expected_y10))
    {
        std::cerr << "p_z basis value on the z-axis is incorrect" << std::endl;
        return 1;
    }
    if (!NearlyEqual(px_axis, 2.0 * expected_y11))
    {
        std::cerr << "p_x basis value on the x-axis is incorrect" << std::endl;
        return 1;
    }
    if (!NearlyEqual(outside_cutoff, 0.0))
    {
        std::cerr << "Basis function should vanish outside the tabulated radial range" << std::endl;
        return 1;
    }

    return 0;
}
