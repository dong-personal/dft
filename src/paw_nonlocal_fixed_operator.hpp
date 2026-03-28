#pragma once

#include "atom.h"
#include "fespace.h"
#include "paw.h"
#include "periodic_kd_tree.hpp"
#include "radial_interpolator.hpp"

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

    void Apply(const mfem::Vector &x_true, mfem::Vector &y_true) const;

    const std::vector<dft::PeriodicNeighbor> &neighbors() const { return m_neighbors; }
    const dft::DenseMatrix &fixed_matrix() const { return m_setup.fixed_nonlocal_correction(); }

  private:
    static Vec3 ToShiftedPoint_(const Vec3 &point, const dft::Structure::LatticeVectors &lattice,
                                const std::array<int, 3> &periodic_image);
    static Vec3 Subtract_(const Vec3 &a, const Vec3 &b);
    static Vec3 FractionalToCartesianShift_(const dft::Structure::LatticeVectors &lattice, const std::array<int, 3> &periodic_image);
    static double Norm_(const Vec3 &v);

    void BuildProjectorTable_();

  private:
    const DFTGLLHexSpace &m_space;
    const dft::PAWSetup &m_setup;
    Vec3 m_atom_center{};
    dft::PeriodicGridPointLocator m_locator;
    std::vector<dft::PeriodicNeighbor> m_neighbors;
    std::vector<dft::RadialInterpolator> m_projector_interpolators;
    std::vector<std::vector<double>> m_projector_values;
};
