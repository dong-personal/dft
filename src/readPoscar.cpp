#include "Structure.h"

#include "constants.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace dft
{
// Declared in Structure.h, implemented here to avoid circular dependency
Structure readPoscar(const std::string &filename)
{
    std::ifstream in(filename);
    if (!in)
        throw std::runtime_error("Failed to open POSCAR: " + filename);

    std::string line;

    // ---------- 1. 注释行 ----------
    std::getline(in, line); // comment, ignore

    // ---------- 2. scale ----------
    std::getline(in, line);
    double scale = std::stod(line);

    // ---------- 3. 晶胞 ----------
    Structure::Mat3 lattice;
    for (int i = 0; i < 3; ++i)
    {
        std::getline(in, line);
        std::istringstream iss(line);
        iss >> lattice[i][0] >> lattice[i][1] >> lattice[i][2];
        for (int j = 0; j < 3; ++j)
            lattice[i][j] *= scale * ANG_TO_BOHR; // 转为 Bohr
    }

    // ---------- 4. 元素符号 ----------
    std::getline(in, line);
    std::istringstream elem_stream(line);
    std::vector<std::string> elements;
    std::string elem;
    while (elem_stream >> elem)
        elements.push_back(elem);

    // ---------- 5. 原子个数 ----------
    std::getline(in, line);
    std::istringstream count_stream(line);
    std::vector<int> counts;
    int c;
    while (count_stream >> c)
        counts.push_back(c);

    if (elements.size() != counts.size())
        throw std::runtime_error("POSCAR: element/count size mismatch");

    // ---------- 6. 坐标类型 ----------
    std::getline(in, line);
    bool fractional = (line[0] == 'D' || line[0] == 'd');

    Structure structure(lattice, fractional ? Structure::CoordType::Fractional : Structure::CoordType::Cartesian);

    // ---------- 7. 原子坐标 ----------
    for (std::size_t i = 0; i < elements.size(); ++i)
    {
        for (int n = 0; n < counts[i]; ++n)
        {
            std::getline(in, line);
            std::istringstream pos_stream(line);

            Atom::Vec3 pos;
            pos_stream >> pos[0] >> pos[1] >> pos[2];

            if (!fractional)
            {
                // Cartesian -> Bohr
                pos[0] *= ANG_TO_BOHR;
                pos[1] *= ANG_TO_BOHR;
                pos[2] *= ANG_TO_BOHR;
            }

            structure.add_atom(Atom(elements[i], pos));
        }
    }

    return structure;
}

} // namespace dft
