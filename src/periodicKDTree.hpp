#pragma once

#include "FEspace.h"
#include "Structure.h"
#include "nanoflann.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace dft
{

struct PeriodicNeighbor
{
    std::size_t index{0};
    double distance_squared{0.0};
    std::array<int, 3> periodic_image{0, 0, 0};
};

class PeriodicKDTree3D
{
  public:
    using Vec3 = std::array<double, 3>;
    using ImageDepth = std::array<int, 3>;

    PeriodicKDTree3D(const std::vector<Vec3> &points, const Structure::Mat3 &lattice,
                     ImageDepth periodic_images = {1, 1, 1});

    std::vector<PeriodicNeighbor> RadiusSearch(const Vec3 &query_point, double radius) const;

    const std::vector<Vec3> &points() const { return points_; }

  private:
    struct BasePoint
    {
        Vec3 position{};
        std::size_t original_index{0};
    };

    struct BasePointCloud
    {
        std::vector<BasePoint> points;

        inline std::size_t kdtree_get_point_count() const { return points.size(); }
        inline double kdtree_get_pt(const std::size_t idx, const std::size_t dim) const
        {
            return points[idx].position[dim];
        }
        template <class BBOX> bool kdtree_get_bbox(BBOX &) const { return false; }
    };

    struct LatticeMatrices
    {
        double cart_from_frac[3][3]{};
        double frac_from_cart[3][3]{};
    };

    static LatticeMatrices BuildLatticeMatrices_(const Structure::Mat3 &lattice);
    static Vec3 WrapToUnitCell_(const Vec3 &point, const LatticeMatrices &matrices, const ImageDepth &periodic_images);
    static Vec3 Multiply_(const double matrix[3][3], const Vec3 &vector);

    void BuildImageShifts_();
    void BuildBasePointCloud_();

  private:
    using KDTreeIndex = nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<double, BasePointCloud>, BasePointCloud, 3, std::size_t>;

    std::vector<Vec3> points_;
    Structure::Mat3 lattice_;
    ImageDepth periodic_images_;
    LatticeMatrices matrices_;
    std::vector<std::array<int, 3>> image_shifts_;
    BasePointCloud base_point_cloud_;
    KDTreeIndex index_;
};

class PeriodicGridPointLocator
{
  public:
    using Vec3 = PeriodicKDTree3D::Vec3;
    using ImageDepth = PeriodicKDTree3D::ImageDepth;

    explicit PeriodicGridPointLocator(const DFTGLLHexSpace &space, ImageDepth periodic_images = {1, 1, 1});

    std::vector<PeriodicNeighbor> RadiusSearch(const Vec3 &query_point, double radius) const;
    const std::vector<Vec3> &points() const { return tree_.points(); }

  private:
    PeriodicKDTree3D tree_;
};

} // namespace dft
