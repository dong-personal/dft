#pragma once
#ifndef DFT_STRUCTURE_H
#define DFT_STRUCTURE_H

#include "atom.h"

#include <algorithm>
#include <vector>

namespace dft
{

class Structure
{
  public:
    using AtomicPosition = Atom::AtomicPosition;
    using LatticeVectors = NDArray<double, 2, int, 9>;

    enum class CoordType
    {
        Cartesian,
        Fractional
    };

    Structure()
        : m_lattice(make_identity_lattice_())
    {
    }

    Structure(const LatticeVectors &lattice, CoordType coord_type = CoordType::Fractional)
        : m_lattice(lattice), m_coord_type(coord_type)
    {
    }

    const LatticeVectors &lattice() const { return m_lattice; }
    void set_lattice(const LatticeVectors &lattice) { m_lattice = lattice; }

    CoordType coord_type() const { return m_coord_type; }
    void set_coord_type(CoordType coord_type) { m_coord_type = coord_type; }

    std::size_t num_atoms() const
    {
        std::size_t total_atoms = 0;
        for (const auto &atom : m_atoms)
        {
            total_atoms += atom.num_positions();
        }
        return total_atoms;
    }

    std::size_t num_species() const { return m_atoms.size(); }

    const std::vector<Atom> &atoms() const { return m_atoms; }
    std::vector<Atom> &atoms() { return m_atoms; }

    void add_atom(const Atom &atom)
    {
        auto atom_it = find_atom_by_symbol_(atom.symbol());
        if (atom_it == m_atoms.end())
        {
            m_atoms.push_back(atom);
            return;
        }

        *atom_it = merge_atoms_(*atom_it, atom);
    }

    void clear_atoms() { m_atoms.clear(); }

  private:
    static LatticeVectors make_identity_lattice_()
    {
        constexpr int kSpatialDimension = 3;
        LatticeVectors lattice(typename LatticeVectors::ShapeType{3, 3});
        for (int row = 0; row < kSpatialDimension; ++row)
        {
            for (int col = 0; col < kSpatialDimension; ++col)
            {
                lattice[row, col] = (row == col) ? 1.0 : 0.0;
            }
        }
        return lattice;
    }

    std::vector<Atom>::iterator find_atom_by_symbol_(const std::string &symbol)
    {
        return std::find_if(m_atoms.begin(), m_atoms.end(),
                            [&](const Atom &existing_atom) { return existing_atom.symbol() == symbol; });
    }

    static Atom merge_atoms_(const Atom &lhs, const Atom &rhs)
    {
        Atom merged_atom(lhs.symbol(), lhs.num_positions() + rhs.num_positions(), lhs.charge());

        for (std::size_t position_index = 0; position_index < lhs.num_positions(); ++position_index)
        {
            merged_atom.set_position(position_index, lhs.position(position_index));
        }

        for (std::size_t position_index = 0; position_index < rhs.num_positions(); ++position_index)
        {
            merged_atom.set_position(lhs.num_positions() + position_index, rhs.position(position_index));
        }

        return merged_atom;
    }

    LatticeVectors m_lattice;
    CoordType m_coord_type{CoordType::Fractional};
    std::vector<Atom> m_atoms;
};

Structure read_poscar(const std::string &filename);
} // namespace dft

#endif // DFT_STRUCTURE_H
