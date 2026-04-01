#pragma once
#ifndef DFT_PKDTREE_HPP
#define DFT_PKDTREE_HPP
#include "fespace.h"
#include "nanoflann.hpp"
#include "structure.h"

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
    using Point3 = DFTMesh::SpatialPoint;
    using ImageDepth = std::array<int, 3>;
    using Matrix3 = Structure::LatticeVectors;

    PeriodicKDTree3D(const mfem::DenseMatrix &point_coordinates, const Structure::LatticeVectors &lattice,
                     ImageDepth periodic_images = {1, 1, 1});
    PeriodicKDTree3D(const std::vector<Point3> &points, const Structure::LatticeVectors &lattice,
                     ImageDepth periodic_images = {1, 1, 1});

    std::vector<PeriodicNeighbor> radius_search(const Point3 &query_point, double radius) const;
    std::vector<PeriodicNeighbor> RadiusSearch(const Point3 &query_point, double radius) const
    {
        return radius_search(query_point, radius);
    }

    Point3 point(std::size_t point_index) const;
    std::size_t point_count() const { return static_cast<std::size_t>(m_points.Width()); }
    const mfem::DenseMatrix &point_coordinates() const { return m_points; }

  private:
    struct BasePointCloud
    {
        mfem::DenseMatrix wrapped_points;
        mfem::Array<int> original_indices;

        BasePointCloud() : wrapped_points(3, 0), original_indices() {}

        std::size_t kdtree_get_point_count() const { return static_cast<std::size_t>(wrapped_points.Width()); }

        double kdtree_get_pt(const std::size_t idx, const std::size_t dim) const
        {
            return wrapped_points(static_cast<int>(dim), static_cast<int>(idx));
        }

        template <class BBOX>
        bool kdtree_get_bbox(BBOX &) const
        {
            return false;
        }
    };

    struct LatticeMatrices
    {
        Matrix3 cartesian_from_fractional{typename Matrix3::ShapeType{3, 3}};
        Matrix3 fractional_from_cartesian{typename Matrix3::ShapeType{3, 3}};
    };

    static mfem::DenseMatrix make_point_matrix(const std::vector<Point3> &points);
    static LatticeMatrices build_lattice_matrices(const Structure::LatticeVectors &lattice);
    static Point3 wrap_to_unit_cell(const Point3 &point, const LatticeMatrices &matrices,
                                    const ImageDepth &periodic_images);
    static Point3 multiply(const Matrix3 &matrix, const Point3 &vector);

    void validate_inputs() const;
    void build_image_shifts();
    void build_base_point_cloud();

  private:
    using KDTreeIndex = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, BasePointCloud>,
                                                            BasePointCloud, 3, std::size_t>;

    mfem::DenseMatrix m_points;
    Structure::LatticeVectors m_lattice;
    ImageDepth m_periodic_images;
    LatticeMatrices m_lattice_matrices;
    std::vector<std::array<int, 3>> m_image_shifts;
    BasePointCloud m_base_point_cloud;
    KDTreeIndex m_index;
};

class PeriodicGridPointLocator
{
  public:
    using Point3 = PeriodicKDTree3D::Point3;
    using ImageDepth = PeriodicKDTree3D::ImageDepth;

    explicit PeriodicGridPointLocator(const DFTGLLHexSpace &space, ImageDepth periodic_images = {1, 1, 1});

    std::vector<PeriodicNeighbor> radius_search(const Point3 &query_point, double radius) const;
    std::vector<PeriodicNeighbor> RadiusSearch(const Point3 &query_point, double radius) const
    {
        return radius_search(query_point, radius);
    }

    Point3 point(std::size_t point_index) const { return m_tree.point(point_index); }
    std::size_t point_count() const { return m_tree.point_count(); }
    const mfem::DenseMatrix &point_coordinates() const { return m_tree.point_coordinates(); }

  private:
    PeriodicKDTree3D m_tree;
};

} // namespace dft

#endif // DFT_PKDTREE_HPP
