#include "fespace.h"
#include "periodic_kd_tree.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

namespace
{

using Vec3 = dft::PeriodicKDTree3D::Vec3;
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
        throw std::runtime_error("singular matrix");
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

Vec3 Multiply(const double matrix[3][3], const Vec3 &vector)
{
    return {
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1] + matrix[0][2] * vector[2],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1] + matrix[1][2] * vector[2],
        matrix[2][0] * vector[0] + matrix[2][1] * vector[1] + matrix[2][2] * vector[2],
    };
}

double DistanceSquared(const Vec3 &x, const Vec3 &y)
{
    const double dx = x[0] - y[0];
    const double dy = x[1] - y[1];
    const double dz = x[2] - y[2];
    return dx * dx + dy * dy + dz * dz;
}

std::vector<NeighborKey> BruteForceNeighborsByImage(const std::vector<Vec3> &points, const dft::Structure::LatticeVectors &lattice,
                                                    const Vec3 &query, double radius,
                                                    const ImageDepth &periodic_images)
{
    double cart_from_frac[3][3];
    for (int row = 0; row < 3; ++row)
    {
        for (int col = 0; col < 3; ++col)
        {
            cart_from_frac[row][col] = lattice[col, row];
        }
    }

    double frac_from_cart[3][3];
    Invert3x3(cart_from_frac, frac_from_cart);

    Vec3 wrapped_query = Multiply(frac_from_cart, query);
    for (int dim = 0; dim < 3; ++dim)
    {
        if (periodic_images[dim] > 0)
        {
            wrapped_query[dim] -= std::floor(wrapped_query[dim]);
        }
    }
    wrapped_query = Multiply(cart_from_frac, wrapped_query);

    const double radius_squared = radius * radius;
    std::vector<NeighborKey> neighbors;
    for (std::size_t i = 0; i < points.size(); ++i)
    {
        Vec3 wrapped_point = Multiply(frac_from_cart, points[i]);
        for (int dim = 0; dim < 3; ++dim)
        {
            if (periodic_images[dim] > 0)
            {
                wrapped_point[dim] -= std::floor(wrapped_point[dim]);
            }
        }
        wrapped_point = Multiply(cart_from_frac, wrapped_point);

        for (int sx = -periodic_images[0]; sx <= periodic_images[0]; ++sx)
        {
            for (int sy = -periodic_images[1]; sy <= periodic_images[1]; ++sy)
            {
                for (int sz = -periodic_images[2]; sz <= periodic_images[2]; ++sz)
                {
                    const Vec3 image_shift = Multiply(cart_from_frac, {static_cast<double>(sx), static_cast<double>(sy),
                                                                       static_cast<double>(sz)});
                    const Vec3 image_point{
                        wrapped_point[0] + image_shift[0],
                        wrapped_point[1] + image_shift[1],
                        wrapped_point[2] + image_shift[2],
                    };
                    if (DistanceSquared(wrapped_query, image_point) <= radius_squared)
                    {
                        neighbors.push_back({i, {sx, sy, sz}});
                    }
                }
            }
        }
    }

    std::sort(neighbors.begin(), neighbors.end(),
              [](const NeighborKey &lhs, const NeighborKey &rhs)
              {
                  if (lhs.periodic_image != rhs.periodic_image)
                  {
                      return lhs.periodic_image < rhs.periodic_image;
                  }
                  return lhs.index < rhs.index;
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
    const auto &points = locator.points();

    if (points.empty())
    {
        std::cerr << "No grid points were extracted from the FE space" << std::endl;
        return 1;
    }

    const Vec3 query = points.front();
    const double radius = 0.8;

    const auto kd_neighbors = locator.RadiusSearch(query, radius);
    auto brute_neighbors = BruteForceNeighborsByImage(points, structure->lattice(), query, radius, periodic_images);

    std::vector<NeighborKey> kd_by_image;
    kd_by_image.reserve(kd_neighbors.size());
    for (const auto &neighbor : kd_neighbors)
    {
        kd_by_image.push_back({neighbor.index, neighbor.periodic_image});
    }

    std::sort(kd_by_image.begin(), kd_by_image.end(),
              [](const NeighborKey &lhs, const NeighborKey &rhs)
              {
                  if (lhs.periodic_image != rhs.periodic_image)
                  {
                      return lhs.periodic_image < rhs.periodic_image;
                  }
                  return lhs.index < rhs.index;
              });

    if (kd_by_image != brute_neighbors)
    {
        std::cerr << "Periodic KD-tree image-aware neighbor list does not match brute-force search" << std::endl;
        return 1;
    }

    const auto &lat = structure->lattice();
    const Vec3 shifted_query{
        query[0] + lat[0, 0],
        query[1] + lat[0, 1],
        query[2] + lat[0, 2],
    };
    const auto shifted_neighbors = locator.RadiusSearch(shifted_query, radius);
    std::vector<NeighborKey> shifted_by_image;
    shifted_by_image.reserve(shifted_neighbors.size());
    for (const auto &neighbor : shifted_neighbors)
    {
        shifted_by_image.push_back({neighbor.index, neighbor.periodic_image});
    }
    std::sort(shifted_by_image.begin(), shifted_by_image.end(),
              [](const NeighborKey &lhs, const NeighborKey &rhs)
              {
                  if (lhs.periodic_image != rhs.periodic_image)
                  {
                      return lhs.periodic_image < rhs.periodic_image;
                  }
                  return lhs.index < rhs.index;
              });

    if (shifted_by_image != kd_by_image)
    {
        std::cerr << "Image-aware periodic KD-tree search should be invariant under lattice translations"
                  << std::endl;
        return 1;
    }

    const ImageDepth mixed_periodicity{1, 0, 0};
    dft::PeriodicGridPointLocator mixed_locator(fespace, mixed_periodicity);
    const auto mixed_kd_neighbors = mixed_locator.RadiusSearch(query, radius);
    const auto mixed_brute = BruteForceNeighborsByImage(points, structure->lattice(), query, radius, mixed_periodicity);

    std::vector<NeighborKey> mixed_kd_by_image;
    mixed_kd_by_image.reserve(mixed_kd_neighbors.size());
    for (const auto &neighbor : mixed_kd_neighbors)
    {
        mixed_kd_by_image.push_back({neighbor.index, neighbor.periodic_image});
    }
    std::sort(mixed_kd_by_image.begin(), mixed_kd_by_image.end(),
              [](const NeighborKey &lhs, const NeighborKey &rhs)
              {
                  if (lhs.periodic_image != rhs.periodic_image)
                  {
                      return lhs.periodic_image < rhs.periodic_image;
                  }
                  return lhs.index < rhs.index;
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
