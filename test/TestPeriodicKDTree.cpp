#include "fespace.h"
#include "pkdtree.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

namespace
{

using Point3 = dft::PeriodicKDTree3D::Point3;
using Matrix3 = dft::PeriodicKDTree3D::Matrix3;
using ImageDepth = dft::PeriodicKDTree3D::ImageDepth;

struct NeighborKey
{
    std::size_t index{0};
    std::array<int, 3> periodic_image{0, 0, 0};

    bool operator==(const NeighborKey &other) const
    {
        return index == other.index && periodic_image == other.periodic_image;
    }
};

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
    if (std::abs(determinant) < 1e-14)
    {
        throw std::runtime_error("singular matrix");
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

Point3 multiply(const Matrix3 &matrix, const Point3 &vector)
{
    return make_point(matrix[0, 0] * vector[0] + matrix[0, 1] * vector[1] + matrix[0, 2] * vector[2],
                      matrix[1, 0] * vector[0] + matrix[1, 1] * vector[1] + matrix[1, 2] * vector[2],
                      matrix[2, 0] * vector[0] + matrix[2, 1] * vector[1] + matrix[2, 2] * vector[2]);
}

double distance_squared(const Point3 &left, const Point3 &right)
{
    const double dx = left[0] - right[0];
    const double dy = left[1] - right[1];
    const double dz = left[2] - right[2];
    return dx * dx + dy * dy + dz * dz;
}

std::vector<NeighborKey> brute_force_neighbors_by_image(const std::vector<Point3> &points,
                                                        const dft::Structure::LatticeVectors &lattice,
                                                        const Point3 &query,
                                                        double radius,
                                                        const ImageDepth &periodic_images)
{
    Matrix3 cartesian_from_fractional(typename Matrix3::ShapeType{3, 3});
    for (int row = 0; row < 3; ++row)
    {
        for (int col = 0; col < 3; ++col)
        {
            cartesian_from_fractional[row, col] = lattice[col, row];
        }
    }

    Matrix3 fractional_from_cartesian(typename Matrix3::ShapeType{3, 3});
    invert_3x3(cartesian_from_fractional, fractional_from_cartesian);

    Point3 wrapped_query = multiply(fractional_from_cartesian, query);
    for (int dim = 0; dim < 3; ++dim)
    {
        if (periodic_images[dim] > 0)
        {
            wrapped_query[dim] -= std::floor(wrapped_query[dim]);
        }
    }
    wrapped_query = multiply(cartesian_from_fractional, wrapped_query);

    const double radius_squared = radius * radius;
    std::vector<NeighborKey> neighbors;
    for (std::size_t point_index = 0; point_index < points.size(); ++point_index)
    {
        Point3 wrapped_point = multiply(fractional_from_cartesian, points[point_index]);
        for (int dim = 0; dim < 3; ++dim)
        {
            if (periodic_images[dim] > 0)
            {
                wrapped_point[dim] -= std::floor(wrapped_point[dim]);
            }
        }
        wrapped_point = multiply(cartesian_from_fractional, wrapped_point);

        for (int shift_x = -periodic_images[0]; shift_x <= periodic_images[0]; ++shift_x)
        {
            for (int shift_y = -periodic_images[1]; shift_y <= periodic_images[1]; ++shift_y)
            {
                for (int shift_z = -periodic_images[2]; shift_z <= periodic_images[2]; ++shift_z)
                {
                    const Point3 image_shift = multiply(cartesian_from_fractional,
                                                        make_point(static_cast<double>(shift_x),
                                                                   static_cast<double>(shift_y),
                                                                   static_cast<double>(shift_z)));
                    const Point3 image_point = make_point(wrapped_point[0] + image_shift[0],
                                                          wrapped_point[1] + image_shift[1],
                                                          wrapped_point[2] + image_shift[2]);
                    if (distance_squared(wrapped_query, image_point) <= radius_squared)
                    {
                        neighbors.push_back({point_index, {shift_x, shift_y, shift_z}});
                    }
                }
            }
        }
    }

    std::sort(neighbors.begin(), neighbors.end(),
              [](const NeighborKey &left, const NeighborKey &right)
              {
                  if (left.periodic_image != right.periodic_image)
                  {
                      return left.periodic_image < right.periodic_image;
                  }
                  return left.index < right.index;
              });

    return neighbors;
}

} // namespace

int main()
{
    const std::string poscar_path = std::string(DFT_SOURCE_DIR) + "/data/C.POSCAR";

    auto structure = std::make_shared<dft::Structure>(dft::read_poscar(poscar_path));
    auto dft_mesh = std::make_shared<dft::DFTMesh>();
    dft_mesh->set_structure(structure);
    dft_mesh->init_periodic_cell_from_lattice(structure->lattice(), 4, 4, 4);

    DFTGLLHexSpace fespace(dft_mesh, 1);
    const ImageDepth periodic_images{1, 1, 1};
    dft::PeriodicGridPointLocator locator(fespace, periodic_images);

    std::vector<Point3> points;
    points.reserve(locator.point_count());
    for (std::size_t point_index = 0; point_index < locator.point_count(); ++point_index)
    {
        points.push_back(locator.point(point_index));
    }

    if (points.empty())
    {
        std::cerr << "No grid points were extracted from the FE space" << std::endl;
        return 1;
    }

    const Point3 query = points.front();
    const double radius = 0.8;

    const auto kd_neighbors = locator.radius_search(query, radius);
    const auto brute_neighbors =
        brute_force_neighbors_by_image(points, structure->lattice(), query, radius, periodic_images);

    std::vector<NeighborKey> kd_by_image;
    kd_by_image.reserve(kd_neighbors.size());
    for (const auto &neighbor : kd_neighbors)
    {
        kd_by_image.push_back({neighbor.index, neighbor.periodic_image});
    }

    std::sort(kd_by_image.begin(), kd_by_image.end(),
              [](const NeighborKey &left, const NeighborKey &right)
              {
                  if (left.periodic_image != right.periodic_image)
                  {
                      return left.periodic_image < right.periodic_image;
                  }
                  return left.index < right.index;
              });

    if (kd_by_image != brute_neighbors)
    {
        std::cerr << "Periodic KD-tree image-aware neighbor list does not match brute-force search" << std::endl;
        return 1;
    }

    const auto &lattice = structure->lattice();
    const Point3 shifted_query = make_point(query[0] + lattice[0, 0], query[1] + lattice[0, 1], query[2] + lattice[0, 2]);
    const auto shifted_neighbors = locator.radius_search(shifted_query, radius);

    std::vector<NeighborKey> shifted_by_image;
    shifted_by_image.reserve(shifted_neighbors.size());
    for (const auto &neighbor : shifted_neighbors)
    {
        shifted_by_image.push_back({neighbor.index, neighbor.periodic_image});
    }

    std::sort(shifted_by_image.begin(), shifted_by_image.end(),
              [](const NeighborKey &left, const NeighborKey &right)
              {
                  if (left.periodic_image != right.periodic_image)
                  {
                      return left.periodic_image < right.periodic_image;
                  }
                  return left.index < right.index;
              });

    if (shifted_by_image != kd_by_image)
    {
        std::cerr << "Image-aware periodic KD-tree search should be invariant under lattice translations"
                  << std::endl;
        return 1;
    }

    const ImageDepth mixed_periodicity{1, 0, 0};
    dft::PeriodicGridPointLocator mixed_locator(fespace, mixed_periodicity);
    const auto mixed_kd_neighbors = mixed_locator.radius_search(query, radius);
    const auto mixed_brute =
        brute_force_neighbors_by_image(points, structure->lattice(), query, radius, mixed_periodicity);

    std::vector<NeighborKey> mixed_kd_by_image;
    mixed_kd_by_image.reserve(mixed_kd_neighbors.size());
    for (const auto &neighbor : mixed_kd_neighbors)
    {
        mixed_kd_by_image.push_back({neighbor.index, neighbor.periodic_image});
    }

    std::sort(mixed_kd_by_image.begin(), mixed_kd_by_image.end(),
              [](const NeighborKey &left, const NeighborKey &right)
              {
                  if (left.periodic_image != right.periodic_image)
                  {
                      return left.periodic_image < right.periodic_image;
                  }
                  return left.index < right.index;
              });

    if (mixed_kd_by_image != mixed_brute)
    {
        std::cerr << "Configurable periodic-image depths do not match brute-force search" << std::endl;
        return 1;
    }

    std::cout << "Grid points: " << points.size() << '\n';
    std::cout << "Image-resolved hits within radius " << radius << ": " << kd_by_image.size() << '\n';
    std::cout << "Mixed-periodicity hits: " << mixed_kd_by_image.size() << '\n';
    std::cout << "First neighbor index: " << kd_neighbors.front().index << '\n';

    return 0;
}
