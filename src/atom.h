#pragma once
#ifndef DFT_ATOM_H
#define DFT_ATOM_H

#include "ndarray/ndarray.hpp"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>

namespace dft
{

class Atom
{
  public:
    using AtomicPosition = NDArray<double, 1, int, 3>;
    using PositionList = NDArray<double, 2, int, 0>;

    Atom() : m_positions(make_position_list_(0)) {}

    // explicit Atom(std::string symbol, double charge = 0.0)
    //     : m_symbol(std::move(symbol)), m_positions(make_position_list_(0)), m_charge(charge)
    // {
    // }

    Atom(std::string symbol, std::size_t num_positions, double charge)
        : m_symbol(std::move(symbol)), m_positions(make_position_list_(num_positions)), m_charge(charge)
    {
    }

    Atom(std::string symbol, std::size_t num_positions) : Atom(std::move(symbol), num_positions, 0.0) {}

    // Atom(std::string symbol, const AtomicPosition &atomic_position, double charge = 0.0)
    //     : Atom(std::move(symbol), 1, charge)
    // {
    //     set_position(0, atomic_position);
    // }

    Atom(std::string symbol, const PositionList &positions, double charge = 0.0)
        : m_symbol(std::move(symbol)), m_positions(positions), m_charge(charge)
    {
        validate_position_list_(m_positions);
    }

    const std::string &symbol() const { return m_symbol; }
    void set_symbol(const std::string &symbol) { m_symbol = symbol; }

    const PositionList &positions() const { return m_positions; }
    PositionList &positions() { return m_positions; }

    std::size_t num_positions() const { return static_cast<std::size_t>(m_positions.shape()[0]); }

    static AtomicPosition read_atomic_position(const PositionList &position_list, std::size_t position_index)
    {
        validate_position_list_(position_list);
        if (position_index >= static_cast<std::size_t>(position_list.shape()[0]))
        {
            throw std::out_of_range("Atom position index out of range");
        }

        AtomicPosition atomic_position(typename AtomicPosition::ShapeType{3});
        for (int dim = 0; dim < 3; ++dim)
        {
            atomic_position[dim] = position_list[static_cast<int>(position_index), dim];
        }
        return atomic_position;
    }

    AtomicPosition position(std::size_t index = 0) const { return read_atomic_position(m_positions, index); }

    void set_position(std::size_t position_index, const AtomicPosition &atomic_position)
    {
        if (position_index >= num_positions())
        {
            throw std::out_of_range("Atom position index out of range");
        }

        for (int dim = 0; dim < 3; ++dim)
        {
            m_positions[static_cast<int>(position_index), dim] = atomic_position[dim];
        }
    }

    double charge() const { return m_charge; }
    void set_charge(double charge) { m_charge = charge; }

  private:
    static PositionList make_position_list_(std::size_t num_positions)
    {
        return PositionList(typename PositionList::ShapeType{static_cast<int>(num_positions), 3});
    }

    static void validate_position_list_(const PositionList &position_list)
    {
        const auto shape = position_list.shape();
        if (shape[1] != 3)
        {
            throw std::invalid_argument("Atom position list must have shape (n, 3)");
        }
    }

  private:
    std::string m_symbol;
    PositionList m_positions;
    double m_charge{0.0};
};

} // namespace dft

#endif // DFT_ATOM_H
