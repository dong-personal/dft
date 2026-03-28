#include "fespace.h"
#include "paw.h"
#include "paw_nonlocal_fixed_operator.hpp"

#include <iostream>
#include <memory>

int main()
{
    const std::string poscar_path = std::string(DFT_SOURCE_DIR) + "/data/C.POSCAR";
    const std::string paw_path = std::string(DFT_SOURCE_DIR) + "/data/C.GGA_PBE-JTH.xml";

    auto structure = std::make_shared<dft::Structure>(dft::read_poscar(poscar_path));
    auto dft_mesh = std::make_shared<dft::DFTMesh>();
    dft_mesh->set_structure(structure);
    dft_mesh->init_periodic_cell_from_lattice(structure->lattice(), 4, 4, 4);

    DFTGLLHexSpace fespace(dft_mesh, 1);
    const dft::PAWSetup setup = dft::LoadPAWSetupXML(paw_path);
    const dft::Atom &atom_species = structure->atoms().front();
    const std::size_t position_index = 0;

    PAWNonlocalFixedOperator nonlocal_fixed(fespace, setup, atom_species, position_index);

    mfem::Vector x(fespace.getDOF());
    x = 1.0;
    mfem::Vector y;
    nonlocal_fixed.Apply(x, y);

    std::cout << "True dofs: " << fespace.getDOF() << '\n';
    std::cout << "Projectors: " << setup.projectors_by_state().size() << '\n';
    std::cout << "Fixed matrix size: " << setup.fixed_nonlocal_correction().size() << '\n';
    std::cout << "Neighbor hits: " << nonlocal_fixed.neighbors().size() << '\n';
    std::cout << "||V_nonlocal^fixed * 1||_2 = " << y.Norml2() << '\n';

    if (setup.fixed_nonlocal_correction().size() != setup.projectors_by_state().size())
    {
        std::cerr << "Fixed nonlocal matrix size does not match projector count" << std::endl;
        return 1;
    }

    if (nonlocal_fixed.neighbors().empty())
    {
        std::cerr << "No nearby grid points were found for the PAW nonlocal operator" << std::endl;
        return 1;
    }

    if (y.Size() != fespace.getDOF())
    {
        std::cerr << "Operator output size mismatch" << std::endl;
        return 1;
    }

    return 0;
}
