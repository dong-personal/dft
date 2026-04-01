#include "paw/ham_correction.hpp"
#include "paw/paw_setup.hpp"

#include <iostream>
#include <string>

namespace
{

std::size_t magnetic_channel_count(const dft::PAWSetup &setup)
{
    std::size_t count = 0;
    for (const dft::PAWState &state : setup.states())
    {
        count += static_cast<std::size_t>(2 * state.l + 1);
    }

    return count;
}

} // namespace

int main()
{
    const std::string paw_path = std::string(DFT_SOURCE_DIR) + "/data/C.GGA_PBE-JTH.xml";

    dft::PAWSetupRegistry registry;
    registry.add(dft::load_paw_setup_xml(paw_path));

    const dft::PAWSetup &stored = registry.get("C");
    const dft::HamCorrection correction = dft::build_ham_correction(stored);

    std::cout << "Stored PAW setup for " << stored.symbol() << '\n';
    std::cout << "Valence charge: " << stored.valence_charge() << '\n';
    std::cout << "Cutoff radius: " << stored.cutoff_radius() << '\n';
    std::cout << "Named metadata blocks: " << stored.metadata_blocks().size() << '\n';
    std::cout << "Named radial grids: " << stored.radial_grids().size() << '\n';
    std::cout << "Named radial functions: " << stored.named_radial_functions().size() << '\n';
    std::cout << "Valence states: " << stored.states().size() << '\n';
    std::cout << "Channels: " << stored.num_channels() << '\n';
    std::cout << "Projectors: " << stored.num_projectors() << '\n';
    std::cout << "Local potential samples: " << stored.local_potential().size() << '\n';
    std::cout << "Core density samples: " << stored.core_density().size() << '\n';
    std::cout << "Fixed nonlocal matrix size: " << correction.fixed_nonlocal_correction().shape()[0] << '\n';

    if (!registry.has("C"))
    {
        std::cerr << "PAW setup registry lookup failed" << std::endl;
        return 1;
    }

    if (stored.num_projectors() != 4)
    {
        std::cerr << "Unexpected projector count: " << stored.num_projectors() << std::endl;
        return 1;
    }

    if (stored.num_channels() != 2)
    {
        std::cerr << "Unexpected channel count: " << stored.num_channels() << std::endl;
        return 1;
    }

    if (stored.local_potential().size() == 0 || stored.core_density().size() == 0)
    {
        std::cerr << "Parsed PAW radial data is empty" << std::endl;
        return 1;
    }

    if (stored.states().size() != 4)
    {
        std::cerr << "Unexpected valence-state count: " << stored.states().size() << std::endl;
        return 1;
    }

    if (stored.radial_grids().find("log1") == stored.radial_grids().end())
    {
        std::cerr << "Missing named radial grid log1" << std::endl;
        return 1;
    }

    if (stored.named_radial_functions().find("blochl_local_ionic_potential") == stored.named_radial_functions().end())
    {
        std::cerr << "Missing blochl_local_ionic_potential dataset" << std::endl;
        return 1;
    }

    if (static_cast<std::size_t>(correction.kinetic_energy_differences().shape()[0]) != magnetic_channel_count(stored))
    {
        std::cerr << "Unexpected kinetic-energy-differences matrix size" << std::endl;
        return 1;
    }

    if (static_cast<std::size_t>(correction.static_coulomb_correction().shape()[0]) != magnetic_channel_count(stored))
    {
        std::cerr << "Unexpected static Coulomb correction matrix size" << std::endl;
        return 1;
    }

    if (static_cast<std::size_t>(correction.fixed_nonlocal_correction().shape()[0]) != magnetic_channel_count(stored))
    {
        std::cerr << "Unexpected fixed nonlocal correction matrix size" << std::endl;
        return 1;
    }

    return 0;
}
