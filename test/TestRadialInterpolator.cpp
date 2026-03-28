#include "paw.h"
#include "radial_interpolator.hpp"

#include <cmath>
#include <iostream>
#include <string>

int main()
{
    const std::string paw_path = std::string(DFT_SOURCE_DIR) + "/data/C.GGA_PBE-JTH.xml";
    const dft::PAWSetup setup = dft::LoadPAWSetupXML(paw_path);

    const auto radial_it = setup.named_radial_functions().find("zero_potential");
    if (radial_it == setup.named_radial_functions().end())
    {
        std::cerr << "Missing zero_potential radial function" << std::endl;
        return 1;
    }

    const dft::RadialFunction &rf = radial_it->second;
    dft::RadialInterpolator interpolator(rf);

    const double v0 = interpolator.Evaluate(rf.r.front());
    const double v1 = interpolator.Evaluate(rf.r[1]);
    const double mid_r = 0.5 * (rf.r[0] + rf.r[1]);
    const double vmid = interpolator.Evaluate(mid_r);
    const double expected_mid = 0.5 * (rf.values[0] + rf.values[1]);

    std::cout << "Zero potential grid points: " << rf.size() << '\n';
    std::cout << "V(r0) = " << v0 << '\n';
    std::cout << "V(r1) = " << v1 << '\n';
    std::cout << "V(mid) = " << vmid << '\n';
    std::cout << "Expected linear midpoint = " << expected_mid << '\n';

    if (std::abs(v0 - rf.values.front()) > 1e-14)
    {
        std::cerr << "Interpolator does not reproduce the first knot" << std::endl;
        return 1;
    }

    if (std::abs(v1 - rf.values[1]) > 1e-14)
    {
        std::cerr << "Interpolator does not reproduce the second knot" << std::endl;
        return 1;
    }

    if (std::abs(vmid - expected_mid) > 1e-12)
    {
        std::cerr << "Linear interpolation midpoint check failed" << std::endl;
        return 1;
    }

    return 0;
}
