#pragma once
#ifndef DFT_DFTMESH_H
#define DFT_DFTMESH_H

#include "structure.h"
#include "mfem.hpp"

#include <memory>
#include <string>
#include <vector>

namespace dft
{

class DFTMesh
{
  public:
    using SpatialPoint = NDArray<double, 1, int, 3>;

    DFTMesh() = default;
    DFTMesh(const DFTMesh &) = delete;
    DFTMesh &operator=(const DFTMesh &) = delete;
    DFTMesh(DFTMesh &&) = delete;
    DFTMesh &operator=(DFTMesh &&) = delete;

    void set_structure(const Structure *structure);
    void set_structure(const Structure &structure);
    void set_structure(std::shared_ptr<const Structure> structure);

    mfem::Mesh &mesh();
    const mfem::Mesh &mesh() const;

    const Structure *structure() const { return m_structure; }
    const Structure::LatticeVectors &lattice() const { return m_lattice_vectors; }

    void save(const std::string &file) const;
    void save_vtu(const std::string &file) const;
    void saveVTU(const std::string &file) const { save_vtu(file); }

    void init_periodic_cell_from_lattice(const Structure::LatticeVectors &lattice_vectors, int nx, int ny, int nz,
                                         mfem::Element::Type element_type = mfem::Element::HEXAHEDRON,
                                         bool use_sfc_ordering = true);

    void init_periodic_cell_from_structure(int nx, int ny, int nz,
                                           mfem::Element::Type element_type = mfem::Element::HEXAHEDRON,
                                           bool use_sfc_ordering = true);

    void refine_near_atoms(double radius, int levels, double shrink = 1.0);

  private:
    static constexpr int kSpatialDimension = 3;

    void ensure_mesh_() const;
    void ensure_base_mesh_() const;
    void ensure_structure_() const;
    void rebuild_periodic_from_base_();

    static std::vector<mfem::Vector>
    make_translations_from_lattice_(const Structure::LatticeVectors &lattice_vectors);
    static mfem::Vector fractional_to_cartesian_(const Structure::LatticeVectors &lattice_vectors,
                                                 const Structure::AtomicPosition &fractional_position);
    static mfem::Vector element_center_from_transform_(mfem::Mesh &mesh, int element_index);
    static double dist2_(const mfem::Vector &lhs, const mfem::Vector &rhs);

    std::vector<mfem::Vector> atom_positions_cartesian_() const;
    double periodic_dist2_(const mfem::Vector &point, const mfem::Vector &atom_position) const;

  private:
    std::unique_ptr<mfem::Mesh> m_mesh;
    std::unique_ptr<mfem::Mesh> m_base_mesh;

    std::shared_ptr<const Structure> m_structure_owner;
    const Structure *m_structure{nullptr};

    Structure::LatticeVectors m_lattice_vectors{typename Structure::LatticeVectors::ShapeType{3, 3}};
    std::vector<mfem::Vector> m_translations;
};

} // namespace dft

#endif // DFT_DFTMESH_H
