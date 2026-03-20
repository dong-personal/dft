#pragma once
#ifndef DFT_STRUCTURE_H
#define DFT_STRUCTURE_H

#include "Atom.h"
#include <array>
#include <vector>

namespace dft
{

class Structure
{
  public:
    using Vec3 = std::array<double, 3>;
    using Mat3 = std::array<Vec3, 3>; // 行向量: a1, a2, a3

    enum class CoordType
    {
        Cartesian,
        Fractional
    };

  public:
    Structure() = default;

    // 构造函数
    Structure(const Mat3 &lattice, CoordType coord_type = CoordType::Fractional)
        : lattice_(lattice), coord_type_(coord_type)
    {
    }

    // ---------- 晶胞 ----------
    const Mat3 &lattice() const { return lattice_; }
    void set_lattice(const Mat3 &lat) { lattice_ = lat; }

    // ---------- 坐标类型 ----------
    CoordType coord_type() const { return coord_type_; }
    void set_coord_type(CoordType t) { coord_type_ = t; }

    // ---------- 原子 ----------
    std::size_t num_atoms() const { return atoms_.size(); }

    const std::vector<Atom> &atoms() const { return atoms_; }
    std::vector<Atom> &atoms() { return atoms_; }

    void add_atom(const Atom &atom) { atoms_.push_back(atom); }

    void clear_atoms() { atoms_.clear(); }

  private:
    Mat3 lattice_{{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};

    CoordType coord_type_{CoordType::Fractional};
    std::vector<Atom> atoms_;
};

// 从 POSCAR 文件读取结构信息 implementated in readPoscar.cpp
Structure readPoscar(const std::string &filename);
} // namespace dft

#endif // DFT_STRUCTURE_H
