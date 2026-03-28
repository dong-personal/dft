#include "paw.h"

#include <iostream>
#include <string>

int main()
{
    const std::string paw_path = std::string(DFT_SOURCE_DIR) + "/data/C.GGA_PBE-JTH.xml";

    dft::PAWSetupRegistry registry;
    registry.Add(dft::LoadPAWSetupXML(paw_path));

    const dft::PAWSetup &stored = registry.Get("C");

    std::cout << "Stored PAW setup for " << stored.symbol() << '\n';
    std::cout << "Valence charge: " << stored.valence_charge() << '\n';
    std::cout << "Cutoff radius: " << stored.cutoff_radius() << '\n';
    std::cout << "Named metadata blocks: " << stored.metadata_blocks().size() << '\n';
    std::cout << "Named radial grids: " << stored.radial_grids().size() << '\n';
    std::cout << "Named radial functions: " << stored.named_radial_functions().size() << '\n';
    std::cout << "Valence states: " << stored.states().size() << '\n';
    std::cout << "Channels: " << stored.NumChannels() << '\n';
    std::cout << "Projectors: " << stored.NumProjectors() << '\n';
    std::cout << "Local potential samples: " << stored.local_potential().size() << '\n';
    std::cout << "Core density samples: " << stored.core_density().size() << '\n';
    std::cout << "Fixed nonlocal matrix size: " << stored.fixed_nonlocal_correction().size() << '\n';

    if (!registry.Has("C"))
    {
        std::cerr << "PAW setup registry lookup failed" << std::endl;
        return 1;
    }

    if (stored.NumProjectors() != 4)
    {
        std::cerr << "Unexpected projector count: " << stored.NumProjectors() << std::endl;
        return 1;
    }

    if (stored.NumChannels() != 2)
    {
        std::cerr << "Unexpected channel count: " << stored.NumChannels() << std::endl;
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

    if (stored.kinetic_energy_differences().size() != stored.states().size())
    {
        std::cerr << "Unexpected kinetic-energy-differences matrix size" << std::endl;
        return 1;
    }

    if (stored.static_coulomb_correction().size() != stored.states().size())
    {
        std::cerr << "Unexpected static Coulomb correction matrix size" << std::endl;
        return 1;
    }

    if (stored.fixed_nonlocal_correction().size() != stored.states().size())
    {
        std::cerr << "Unexpected fixed nonlocal correction matrix size" << std::endl;
        return 1;
    }

    if (stored.channels().front().fixed_nonlocal_correction.size() != stored.channels().front().projectors.size())
    {
        std::cerr << "Channel fixed nonlocal matrix does not match projector count" << std::endl;
        return 1;
    }

    return 0;
}
