#include "Structure.h"
#include <iomanip>
#include <iostream>

int main() {
  try {
    // 1. 读取 POSCAR
    dft::Structure structure = dft::readPoscar("C.POSCAR");

    // 2. 打印晶胞
    std::cout << "Lattice vectors:\n";
    const auto &lat = structure.lattice();
    for (int i = 0; i < 3; ++i) {
      std::cout << std::setw(12) << lat[i][0] << std::setw(12) << lat[i][1]
                << std::setw(12) << lat[i][2] << "\n";
    }

    // 3. 打印坐标类型
    std::cout << "\nCoordinate type: ";
    if (structure.coord_type() == dft::Structure::CoordType::Fractional)
      std::cout << "Fractional\n";
    else
      std::cout << "Cartesian\n";

    // 4. 打印原子信息
    std::cout << "\nAtoms (" << structure.num_atoms() << "):\n";
    for (std::size_t i = 0; i < structure.atoms().size(); ++i) {
      const auto &atom = structure.atoms()[i];
      const auto &p = atom.position();

      std::cout << std::setw(3) << i << "  " << std::setw(2) << atom.symbol()
                << "  " << std::fixed << std::setprecision(6) << std::setw(10)
                << p[0] << std::setw(10) << p[1] << std::setw(10) << p[2]
                << "\n";
    }

    std::cout << "\nPOSCAR read successfully.\n";
  } catch (const std::exception &e) {
    std::cerr << "Error reading POSCAR: " << e.what() << "\n";
    return 1;
  }

  return 0;
}
