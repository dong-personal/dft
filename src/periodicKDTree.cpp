#include "periodicKDTree.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace dft
{

namespace
{

double Determinant3x3(const double a[3][3])
{
    return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
           a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
           a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
}

void Invert3x3(const double a[3][3], double inv[3][3])
{
    const double det = Determinant3x3(a);
    if (std::abs(det) < 1e-14)
    {
        throw std::runtime_error("PeriodicKDTree3D requires a non-singular lattice");
    }

    inv[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]) / det;
    inv[0][1] = (a[0][2] * a[2][1] - a[0][1] * a[2][2]) / det;
    inv[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) / det;

    inv[1][0] = (a[1][2] * a[2][0] - a[1][0] * a[2][2]) / det;
    inv[1][1] = (a[0][0] * a[2][2] - a[0][2] * a[2][0]) / det;
    inv[1][2] = (a[0][2] * a[1][0] - a[0][0] * a[1][2]) / det;

    inv[2][0] = (a[1][0] * a[2][1] - a[1][1] * a[2][0]) / det;
    inv[2][1] = (a[0][1] * a[2][0] - a[0][0] * a[2][1]) / det;
    inv[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) / det;
}

} // namespace

PeriodicKDTree3D::PeriodicKDTree3D(const std::vector<Vec3> &points, const Structure::Mat3 &lattice,
                                   ImageDepth periodic_images)
    : points_(points),
      lattice_(lattice),
      periodic_images_(periodic_images),
      matrices_(BuildLatticeMatrices_(lattice_)),
      image_shifts_(),
      base_point_cloud_(),
      index_(3, base_point_cloud_, nanoflann::KDTreeSingleIndexAdaptorParams(10))
{
    for (const int value : periodic_images_)
    {
        if (value < 0)
        {
            throw std::invalid_argument("PeriodicKDTree3D requires non-negative image depths");
        }
    }
    BuildImageShifts_();
    BuildBasePointCloud_();
    index_.buildIndex();
}

std::vector<PeriodicNeighbor> PeriodicKDTree3D::RadiusSearch(const Vec3 &query_point, double radius) const
{
    if (radius < 0.0)
    {
        throw std::invalid_argument("PeriodicKDTree3D::RadiusSearch requires radius >= 0");
    }

    const Vec3 wrapped_query = WrapToUnitCell_(query_point, matrices_, periodic_images_);
    std::vector<PeriodicNeighbor> neighbors;

    nanoflann::SearchParameters search_params;
    for (const auto &image_shift_int : image_shifts_)
    {
        const Vec3 image_shift{
            static_cast<double>(image_shift_int[0]),
            static_cast<double>(image_shift_int[1]),
            static_cast<double>(image_shift_int[2]),
        };
        const Vec3 cartesian_shift = Multiply_(matrices_.cart_from_frac, image_shift);
        const Vec3 shifted_query{
            wrapped_query[0] - cartesian_shift[0],
            wrapped_query[1] - cartesian_shift[1],
            wrapped_query[2] - cartesian_shift[2],
        };

        std::vector<nanoflann::ResultItem<std::size_t, double>> matches;
        index_.radiusSearch(shifted_query.data(), radius * radius, matches, search_params);
        for (const auto &match : matches)
        {
            const BasePoint &base_point = base_point_cloud_.points[match.first];
            neighbors.push_back({base_point.original_index, match.second, image_shift_int});
        }
    }

    std::sort(neighbors.begin(), neighbors.end(),
              [](const PeriodicNeighbor &lhs, const PeriodicNeighbor &rhs)
              {
                  if (lhs.periodic_image != rhs.periodic_image)
                  {
                      return lhs.periodic_image < rhs.periodic_image;
                  }
                  if (lhs.distance_squared != rhs.distance_squared)
                  {
                      return lhs.distance_squared < rhs.distance_squared;
                  }
                  return lhs.index < rhs.index;
              });

    neighbors.erase(std::unique(neighbors.begin(), neighbors.end(),
                                [](const PeriodicNeighbor &lhs, const PeriodicNeighbor &rhs)
                                {
                                    return lhs.index == rhs.index && lhs.periodic_image == rhs.periodic_image;
                                }),
                    neighbors.end());

    return neighbors;
}

PeriodicKDTree3D::LatticeMatrices PeriodicKDTree3D::BuildLatticeMatrices_(const Structure::Mat3 &lattice)
{
    LatticeMatrices matrices;
    for (int row = 0; row < 3; ++row)
    {
        for (int col = 0; col < 3; ++col)
        {
            matrices.cart_from_frac[row][col] = lattice[col][row];
        }
    }
    Invert3x3(matrices.cart_from_frac, matrices.frac_from_cart);
    return matrices;
}

PeriodicKDTree3D::Vec3 PeriodicKDTree3D::WrapToUnitCell_(const Vec3 &point, const LatticeMatrices &matrices,
                                                         const ImageDepth &periodic_images)
{
    Vec3 fractional = Multiply_(matrices.frac_from_cart, point);
    for (int i = 0; i < 3; ++i)
    {
        if (periodic_images[i] > 0)
        {
            fractional[i] -= std::floor(fractional[i]);
        }
    }
    return Multiply_(matrices.cart_from_frac, fractional);
}

PeriodicKDTree3D::Vec3 PeriodicKDTree3D::Multiply_(const double matrix[3][3], const Vec3 &vector)
{
    return {
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1] + matrix[0][2] * vector[2],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1] + matrix[1][2] * vector[2],
        matrix[2][0] * vector[0] + matrix[2][1] * vector[1] + matrix[2][2] * vector[2],
    };
}

void PeriodicKDTree3D::BuildImageShifts_()
{
    image_shifts_.clear();
    for (int i = -periodic_images_[0]; i <= periodic_images_[0]; ++i)
    {
        for (int j = -periodic_images_[1]; j <= periodic_images_[1]; ++j)
        {
            for (int k = -periodic_images_[2]; k <= periodic_images_[2]; ++k)
            {
                image_shifts_.push_back({i, j, k});
            }
        }
    }
}

void PeriodicKDTree3D::BuildBasePointCloud_()
{
    base_point_cloud_.points.clear();
    base_point_cloud_.points.reserve(points_.size());

    for (std::size_t original_index = 0; original_index < points_.size(); ++original_index)
    {
        const Vec3 wrapped_point = WrapToUnitCell_(points_[original_index], matrices_, periodic_images_);
        base_point_cloud_.points.push_back({wrapped_point, original_index});
    }
}

PeriodicGridPointLocator::PeriodicGridPointLocator(const DFTGLLHexSpace &space, ImageDepth periodic_images)
    : tree_(space.TrueDofCoordinates(), space.MeshSource().lattice(), periodic_images)
{
}

std::vector<PeriodicNeighbor> PeriodicGridPointLocator::RadiusSearch(const Vec3 &query_point, double radius) const
{
    return tree_.RadiusSearch(query_point, radius);
}

} // namespace dft
