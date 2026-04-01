#pragma once

#ifndef DFT_PAW_GRID_MAP_HPP
#define DFT_PAW_GRID_MAP_HPP
#include "fespace.h"
#include "mfem.hpp"
#include "paw/ham_correction.hpp"
#include "paw/paw_setup.hpp"
#include "paw/projector.hpp"
#include "pkdtree.hpp"

#include <cstddef>
#include <vector>

namespace dft
{

class AtomGLLHamCorrectionMatrix
{
  public:
    AtomGLLHamCorrectionMatrix(const DFTGLLHexSpace &space, const PAWSetup &setup, const HamCorrection &ham_correction,
                               const Atom &atom, std::size_t position_index = 0);

    const DenseMatrix &projector_basis_overlap() const { return m_projector_basis_overlap; }
    const DenseMatrix &fixed_nonlocal() const { return m_fixed_nonlocal; }
    const DenseMatrix &matrix() const { return m_matrix; }

  private:
    DenseMatrix m_projector_basis_overlap;
    DenseMatrix m_fixed_nonlocal;
    DenseMatrix m_matrix;
};

struct AtomGridMap
{
    DenseMatrix values;
    std::vector<std::size_t> point_indices;

    std::size_t point_count() const { return point_indices.size(); }
};

AtomGridMap build_atom_grid_map(const HamCorrection &ham_correction, const ProjectorTable &projector_table,
                                const std::vector<PeriodicNeighbor> &neighbors, const mfem::Vector &mass_diagonal);

DenseMatrix build_atom_gll_ham_correction_matrix(const DFTGLLHexSpace &space, const PAWSetup &setup,
                                                 const HamCorrection &ham_correction, const Atom &atom,
                                                 std::size_t position_index = 0);

DenseMatrix build_gll_ham_correction_matrix(const DFTGLLHexSpace &space, const PAWSetupRegistry &setup_registry);

} // namespace dft

#endif // DFT_PAW_GRID_MAP_HPP
