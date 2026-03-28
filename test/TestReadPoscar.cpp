#include "structure.h"
#include <iomanip>
#include <iostream>

int main() {
  try {
    const std::string poscar_path = std::string(DFT_SOURCE_DIR) + "/data/C.POSCAR";
    dft::Structure structure = dft::read_poscar(poscar_path);

    // 2. 打印晶胞
    std::cout << "Lattice vectors:\n";
    const auto &lat = structure.lattice();
    for (int i = 0; i < 3; ++i) {
      std::cout << std::setw(12) << lat[i, 0] << std::setw(12) << lat[i, 1]
                << std::setw(12) << lat[i, 2] << "\n";
    }

    // 3. 打印坐标类型
    std::cout << "\nCoordinate type: ";
    if (structure.coord_type() == dft::Structure::CoordType::Fractional)
      std::cout << "Fractional\n";
    else
      std::cout << "Cartesian\n";

    // 4. 打印原子信息
    std::cout << "\nAtoms (" << structure.num_atoms() << " total, " << structure.num_species() << " species):\n";
    std::size_t atom_index = 0;
    for (const auto &atom : structure.atoms()) {
      for (std::size_t position_index = 0; position_index < atom.num_positions(); ++position_index) {
        const auto p = atom.position(position_index);
        std::cout << std::setw(3) << atom_index << "  " << std::setw(2) << atom.symbol()
                  << "  " << std::fixed << std::setprecision(6) << std::setw(10)
                  << p[0] << std::setw(10) << p[1] << std::setw(10) << p[2]
                  << "\n";
        ++atom_index;
      }
    }

    std::cout << "\nPOSCAR read successfully.\n";
  } catch (const std::exception &e) {
    std::cerr << "Error reading POSCAR: " << e.what() << "\n";
    return 1;
  }

  return 0;
}
