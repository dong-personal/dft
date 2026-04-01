#include "pkdtree.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dft
{

namespace
{

constexpr int kSpatialDimension = 3;
constexpr int kKdTreeLeafSize = 10;
constexpr double kSingularTolerance = 1e-14;

using Point3 = PeriodicKDTree3D::Point3;
using Matrix3 = PeriodicKDTree3D::Matrix3;

Point3 make_point(double x, double y, double z)
{
    Point3 point(typename Point3::ShapeType{3});
    point[0] = x;
    point[1] = y;
    point[2] = z;
    return point;
}

double determinant_3x3(const Matrix3 &matrix)
{
    return matrix[0, 0] * (matrix[1, 1] * matrix[2, 2] - matrix[1, 2] * matrix[2, 1]) -
           matrix[0, 1] * (matrix[1, 0] * matrix[2, 2] - matrix[1, 2] * matrix[2, 0]) +
           matrix[0, 2] * (matrix[1, 0] * matrix[2, 1] - matrix[1, 1] * matrix[2, 0]);
}

void invert_3x3(const Matrix3 &matrix, Matrix3 &inverse)
{
    const double determinant = determinant_3x3(matrix);
    if (std::abs(determinant) < kSingularTolerance)
    {
        throw std::runtime_error("PeriodicKDTree3D requires a non-singular lattice");
    }

    inverse[0, 0] = (matrix[1, 1] * matrix[2, 2] - matrix[1, 2] * matrix[2, 1]) / determinant;
    inverse[0, 1] = (matrix[0, 2] * matrix[2, 1] - matrix[0, 1] * matrix[2, 2]) / determinant;
    inverse[0, 2] = (matrix[0, 1] * matrix[1, 2] - matrix[0, 2] * matrix[1, 1]) / determinant;

    inverse[1, 0] = (matrix[1, 2] * matrix[2, 0] - matrix[1, 0] * matrix[2, 2]) / determinant;
    inverse[1, 1] = (matrix[0, 0] * matrix[2, 2] - matrix[0, 2] * matrix[2, 0]) / determinant;
    inverse[1, 2] = (matrix[0, 2] * matrix[1, 0] - matrix[0, 0] * matrix[1, 2]) / determinant;

    inverse[2, 0] = (matrix[1, 0] * matrix[2, 1] - matrix[1, 1] * matrix[2, 0]) / determinant;
    inverse[2, 1] = (matrix[0, 1] * matrix[2, 0] - matrix[0, 0] * matrix[2, 1]) / determinant;
    inverse[2, 2] = (matrix[0, 0] * matrix[1, 1] - matrix[0, 1] * matrix[1, 0]) / determinant;
}

mfem::DenseMatrix point_matrix_from_spatial_points(const std::vector<dft::DFTMesh::SpatialPoint> &points)
{
    mfem::DenseMatrix point_matrix(kSpatialDimension, static_cast<int>(points.size()));
    for (int point_index = 0; point_index < point_matrix.Width(); ++point_index)
    {
        for (int dimension = 0; dimension < kSpatialDimension; ++dimension)
        {
            point_matrix(dimension, point_index) = points[static_cast<std::size_t>(point_index)][dimension];
        }
    }

    return point_matrix;
}

} // namespace

PeriodicKDTree3D::PeriodicKDTree3D(const mfem::DenseMatrix &point_coordinates, const Structure::LatticeVectors &lattice,
                                   ImageDepth periodic_images)
    : m_points(point_coordinates),
      m_lattice(lattice),
      m_periodic_images(periodic_images),
      m_lattice_matrices(build_lattice_matrices(m_lattice)),
      m_image_shifts(),
      m_base_point_cloud(),
      m_index(kSpatialDimension, m_base_point_cloud, nanoflann::KDTreeSingleIndexAdaptorParams(kKdTreeLeafSize))
{
    validate_inputs();
    build_image_shifts();
    build_base_point_cloud();
    m_index.buildIndex();
}

PeriodicKDTree3D::PeriodicKDTree3D(const std::vector<Point3> &points, const Structure::LatticeVectors &lattice,
                                   ImageDepth periodic_images)
    : PeriodicKDTree3D(make_point_matrix(points), lattice, periodic_images)
{
}

std::vector<PeriodicNeighbor> PeriodicKDTree3D::radius_search(const Point3 &query_point, double radius) const
{
    if (radius < 0.0)
    {
        throw std::invalid_argument("PeriodicKDTree3D::radius_search requires radius >= 0");
    }

    const Point3 wrapped_query = wrap_to_unit_cell(query_point, m_lattice_matrices, m_periodic_images);
    std::vector<PeriodicNeighbor> neighbors;

    nanoflann::SearchParameters search_parameters;
    for (const std::array<int, 3> &image_shift : m_image_shifts)
    {
        const Point3 fractional_shift = make_point(static_cast<double>(image_shift[0]),
                                                   static_cast<double>(image_shift[1]),
                                                   static_cast<double>(image_shift[2]));
        const Point3 cartesian_shift = multiply(m_lattice_matrices.cartesian_from_fractional, fractional_shift);
        const Point3 shifted_query = make_point(wrapped_query[0] - cartesian_shift[0],
                                                wrapped_query[1] - cartesian_shift[1],
                                                wrapped_query[2] - cartesian_shift[2]);

        std::vector<nanoflann::ResultItem<std::size_t, double>> matches;
        m_index.radiusSearch(shifted_query.data_pointer(), radius * radius, matches, search_parameters);
        for (const auto &match : matches)
        {
            neighbors.push_back(
                {static_cast<std::size_t>(m_base_point_cloud.original_indices[static_cast<int>(match.first)]),
                 match.second,
                 image_shift});
        }
    }

    std::sort(neighbors.begin(), neighbors.end(),
              [](const PeriodicNeighbor &left, const PeriodicNeighbor &right)
              {
                  if (left.periodic_image != right.periodic_image)
                  {
                      return left.periodic_image < right.periodic_image;
                  }
                  if (left.distance_squared != right.distance_squared)
                  {
                      return left.distance_squared < right.distance_squared;
                  }
                  return left.index < right.index;
              });

    neighbors.erase(std::unique(neighbors.begin(), neighbors.end(),
                                [](const PeriodicNeighbor &left, const PeriodicNeighbor &right)
                                {
                                    return left.index == right.index && left.periodic_image == right.periodic_image;
                                }),
                    neighbors.end());

    return neighbors;
}

PeriodicKDTree3D::Point3 PeriodicKDTree3D::point(std::size_t point_index) const
{
    if (point_index >= point_count())
    {
        throw std::out_of_range("PeriodicKDTree3D point index out of range");
    }

    return make_point(m_points(0, static_cast<int>(point_index)),
                      m_points(1, static_cast<int>(point_index)),
                      m_points(2, static_cast<int>(point_index)));
}

mfem::DenseMatrix PeriodicKDTree3D::make_point_matrix(const std::vector<Point3> &points)
{
    mfem::DenseMatrix point_matrix(kSpatialDimension, static_cast<int>(points.size()));
    for (int point_index = 0; point_index < point_matrix.Width(); ++point_index)
    {
        for (int dimension = 0; dimension < kSpatialDimension; ++dimension)
        {
            point_matrix(dimension, point_index) = points[static_cast<std::size_t>(point_index)][dimension];
        }
    }

    return point_matrix;
}

PeriodicKDTree3D::LatticeMatrices PeriodicKDTree3D::build_lattice_matrices(const Structure::LatticeVectors &lattice)
{
    LatticeMatrices matrices;
    for (int row = 0; row < kSpatialDimension; ++row)
    {
        for (int column = 0; column < kSpatialDimension; ++column)
        {
            matrices.cartesian_from_fractional[row, column] = lattice[column, row];
        }
    }

    invert_3x3(matrices.cartesian_from_fractional, matrices.fractional_from_cartesian);
    return matrices;
}

PeriodicKDTree3D::Point3 PeriodicKDTree3D::wrap_to_unit_cell(const Point3 &point,
                                                             const LatticeMatrices &matrices,
                                                             const ImageDepth &periodic_images)
{
    Point3 fractional_coordinates = multiply(matrices.fractional_from_cartesian, point);
    for (int dimension = 0; dimension < kSpatialDimension; ++dimension)
    {
        if (periodic_images[dimension] > 0)
        {
            fractional_coordinates[dimension] -= std::floor(fractional_coordinates[dimension]);
        }
    }

    return multiply(matrices.cartesian_from_fractional, fractional_coordinates);
}

PeriodicKDTree3D::Point3 PeriodicKDTree3D::multiply(const Matrix3 &matrix, const Point3 &vector)
{
    return make_point(matrix[0, 0] * vector[0] + matrix[0, 1] * vector[1] + matrix[0, 2] * vector[2],
                      matrix[1, 0] * vector[0] + matrix[1, 1] * vector[1] + matrix[1, 2] * vector[2],
                      matrix[2, 0] * vector[0] + matrix[2, 1] * vector[1] + matrix[2, 2] * vector[2]);
}

void PeriodicKDTree3D::validate_inputs() const
{
    if (m_points.Height() != kSpatialDimension)
    {
        throw std::invalid_argument("PeriodicKDTree3D requires point coordinates stored as a 3 x N matrix");
    }

    for (const int periodic_depth : m_periodic_images)
    {
        if (periodic_depth < 0)
        {
            throw std::invalid_argument("PeriodicKDTree3D requires non-negative image depths");
        }
    }
}

void PeriodicKDTree3D::build_image_shifts()
{
    m_image_shifts.clear();
    for (int shift_x = -m_periodic_images[0]; shift_x <= m_periodic_images[0]; ++shift_x)
    {
        for (int shift_y = -m_periodic_images[1]; shift_y <= m_periodic_images[1]; ++shift_y)
        {
            for (int shift_z = -m_periodic_images[2]; shift_z <= m_periodic_images[2]; ++shift_z)
            {
                m_image_shifts.push_back({shift_x, shift_y, shift_z});
            }
        }
    }
}

void PeriodicKDTree3D::build_base_point_cloud()
{
    m_base_point_cloud.wrapped_points.SetSize(kSpatialDimension, m_points.Width());
    m_base_point_cloud.original_indices.SetSize(m_points.Width());

    for (int point_index = 0; point_index < m_points.Width(); ++point_index)
    {
        const Point3 point_coordinates =
            make_point(m_points(0, point_index), m_points(1, point_index), m_points(2, point_index));
        const Point3 wrapped_point = wrap_to_unit_cell(point_coordinates, m_lattice_matrices, m_periodic_images);
        for (int dimension = 0; dimension < kSpatialDimension; ++dimension)
        {
            m_base_point_cloud.wrapped_points(dimension, point_index) = wrapped_point[dimension];
        }
        m_base_point_cloud.original_indices[point_index] = point_index;
    }
}

PeriodicGridPointLocator::PeriodicGridPointLocator(const DFTGLLHexSpace &space, ImageDepth periodic_images)
    : m_tree(point_matrix_from_spatial_points(space.true_dof_coordinates()), space.mesh_source().lattice(),
             periodic_images)
{
}

std::vector<PeriodicNeighbor> PeriodicGridPointLocator::radius_search(const Point3 &query_point, double radius) const
{
    return m_tree.radius_search(query_point, radius);
}

} // namespace dft
