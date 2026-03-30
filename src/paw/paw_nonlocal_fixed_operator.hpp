#pragma once

#include "atom.h"
#include "fespace.h"
#include "paw/paw.hpp"
#include "paw/paw_basis_evaluator.hpp"
#include "paw/radial_interpolator.hpp"
#include "periodic_kd_tree.hpp"

#include <array>
#include <cstddef>
#include <vector>

class PAWNonlocalFixedOperator
{
  public:
    using Vec3 = std::array<double, 3>;

    PAWNonlocalFixedOperator(const DFTGLLHexSpace &space, const dft::PAWSetup &setup, const dft::Atom &atom,
                             std::size_t position_index = 0,
                             dft::PeriodicKDTree3D::ImageDepth periodic_images = {1, 1, 1});

    void apply(const mfem::Vector &x_true, mfem::Vector &y_true) const;

    const std::vector<dft::PeriodicNeighbor> &neighbors() const { return m_neighbors; }
    const dft::DenseMatrix &fixed_matrix() const { return m_setup.fixed_nonlocal_correction(); }

  private:
    static Vec3 to_shifted_point(const Vec3 &point, const dft::Structure::LatticeVectors &lattice,
                                 const std::array<int, 3> &periodic_image);
    static Vec3 subtract_vectors(const Vec3 &left, const Vec3 &right);
    static Vec3 fractional_to_cartesian_shift(const dft::Structure::LatticeVectors &lattice,
                                              const std::array<int, 3> &periodic_image);
    static double norm(const Vec3 &vector);

    void build_projector_table();

  private:
    const DFTGLLHexSpace &m_space;
    const dft::PAWSetup &m_setup;
    Vec3 m_atom_center{};
    dft::PeriodicGridPointLocator m_locator;
    std::vector<dft::PeriodicNeighbor> m_neighbors;
    std::vector<dft::RadialInterpolator> m_projector_interpolators;
    std::vector<std::vector<double>> m_projector_values;
};
