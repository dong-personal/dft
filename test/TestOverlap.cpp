#include "FEspace.h"

int main(void)
{
    dft::Structure structure = dft::readPoscar("C.POSCAR");
    dft::DFTMesh dft_mesh;
    dft_mesh.set_structure(&structure);
    dft_mesh.init_periodic_cell_from_lattice(structure.lattice(), 10, 10, 10);
    dft_mesh.refine_near_atoms(1.0, 2, 0.5);

    DFTGLLHexSpace fespace(dft_mesh, 3);

    mfem::Vector diag = fespace.MassDiagTrue();
}