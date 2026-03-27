#pragma once

#include "Atom.h"
#include "FEspace.h"
#include "PAW.h"
#include "periodicKDTree.hpp"
#include "radialInterpolator.hpp"

#include <array>
#include <vector>

class PAWNonlocalFixedOperator
{
  public:
    using Vec3 = std::array<double, 3>;

    PAWNonlocalFixedOperator(const DFTGLLHexSpace &space, const dft::PAWSetup &setup, const dft::Atom &atom,
                             dft::PeriodicKDTree3D::ImageDepth periodic_images = {1, 1, 1});

    void Apply(const mfem::Vector &x_true, mfem::Vector &y_true) const;

    const std::vector<dft::PeriodicNeighbor> &neighbors() const { return neighbors_; }
    const dft::DenseMatrix &fixed_matrix() const { return setup_.fixed_nonlocal_correction(); }

  private:
    static Vec3 ToShiftedPoint_(const Vec3 &point, const dft::Structure::Mat3 &lattice,
                                const std::array<int, 3> &periodic_image);
    static Vec3 Subtract_(const Vec3 &a, const Vec3 &b);
    static Vec3 FractionalToCartesianShift_(const dft::Structure::Mat3 &lattice, const std::array<int, 3> &periodic_image);
    static double Norm_(const Vec3 &v);

    void BuildProjectorTable_();

  private:
    const DFTGLLHexSpace &space_;
    const dft::PAWSetup &setup_;
    dft::Atom atom_;
    dft::PeriodicGridPointLocator locator_;
    std::vector<dft::PeriodicNeighbor> neighbors_;
    std::vector<dft::RadialInterpolator> projector_interpolators_;
    std::vector<std::vector<double>> projector_values_;
};
