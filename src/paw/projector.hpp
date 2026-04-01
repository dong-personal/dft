#pragma once
#ifndef DFT_PAW_PROJECTOR_HPP
#define DFT_PAW_PROJECTOR_HPP
#include "atom.h"
#include "fespace.h"
#include "mfem.hpp"
#include "paw/interpolator.hpp"
#include "paw/paw_setup.hpp"
#include "pkdtree.hpp"

#include <vector>

namespace dft
{

using ProjectorValueTable = NDArray<double, 2, int>;

struct ProjectorTable
{
    ProjectorValueTable values;
    std::vector<std::size_t> state_indices;
    std::vector<int> magnetic_quantum_numbers;

    std::size_t channel_count() const { return state_indices.size(); }
};

ProjectorTable build_projector_table(const PAWSetup &setup, const Structure::LatticeVectors &lattice,
                                     const Atom::AtomicPosition &atom_center, const PeriodicGridPointLocator &locator,
                                     const std::vector<PeriodicNeighbor> &neighbors);

DenseMatrix build_projector_basis_overlap_matrix(const DFTGLLHexSpace &space, const PAWSetup &setup, const Atom &atom,
                                                 std::size_t position_index);

} // namespace dft
#endif // DFT_PAW_PROJECTOR_HPP
