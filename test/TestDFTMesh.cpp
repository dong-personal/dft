#include "dftmesh.h"

#include <string>

int main()
{
  const std::string poscar_path = std::string(DFT_SOURCE_DIR) + "/data/C.POSCAR";

  dft::Structure structure = dft::read_poscar(poscar_path);
  dft::DFTMesh dft_mesh;

  dft_mesh.set_structure(&structure);
  dft_mesh.init_periodic_cell_from_lattice(structure.lattice(), 10, 10, 10);
  dft_mesh.refine_near_atoms(1.0, 2, 0.5);
  dft_mesh.save_vtu("Cmesh.vtu");

  return 0;
}
