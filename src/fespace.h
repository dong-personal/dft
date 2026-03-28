#pragma once
#ifndef DFT_FESPACE_H
#define DFT_FESPACE_H

#include "dftmesh.h"
#include "mfem.hpp"

#include <memory>
#include <vector>

class DFTGLLHexSpace
{
  public:
    using SpatialPoint = dft::DFTMesh::SpatialPoint;

    explicit DFTGLLHexSpace(std::shared_ptr<dft::DFTMesh> mesh, int order, int vdim = 1);
    DFTGLLHexSpace(dft::DFTMesh &mesh, int order, int vdim = 1);

    DFTGLLHexSpace(const DFTGLLHexSpace &) = delete;
    DFTGLLHexSpace &operator=(const DFTGLLHexSpace &) = delete;
    DFTGLLHexSpace(DFTGLLHexSpace &&) = delete;
    DFTGLLHexSpace &operator=(DFTGLLHexSpace &&) = delete;

    mfem::Mesh &mesh() { return *m_mesh; }
    const mfem::Mesh &mesh() const { return *m_mesh; }
    mfem::Mesh &Mesh() { return mesh(); }
    const mfem::Mesh &Mesh() const { return mesh(); }

    dft::DFTMesh &mesh_source() { return *m_dft_mesh; }
    const dft::DFTMesh &mesh_source() const { return *m_dft_mesh; }
    dft::DFTMesh &MeshSource() { return mesh_source(); }
    const dft::DFTMesh &MeshSource() const { return mesh_source(); }

    mfem::FiniteElementSpace &space() { return m_fes; }
    const mfem::FiniteElementSpace &space() const { return m_fes; }
    mfem::FiniteElementSpace &Space() { return space(); }
    const mfem::FiniteElementSpace &Space() const { return space(); }

    int order() const { return m_order; }
    int Order() const { return order(); }

    int true_dof_count() const { return m_fes.GetTrueVSize(); }
    int getDOF() const { return true_dof_count(); }

    const mfem::Vector &mass_diag_true() const { return m_mass_diag_true; }
    const mfem::Vector &MassDiagTrue() const { return mass_diag_true(); }

    const mfem::Vector &minv_half_diag_true() const { return m_minv_half_diag_true; }
    const mfem::Vector &MinvHalfDiagTrue() const { return minv_half_diag_true(); }

    std::vector<SpatialPoint> true_dof_coordinates() const;
    std::vector<SpatialPoint> TrueDofCoordinates() const { return true_dof_coordinates(); }

    void apply_mass_diag_true(const mfem::Vector &x_true, mfem::Vector &y_true) const;
    void ApplyMassDiagTrue(const mfem::Vector &x_true, mfem::Vector &y_true) const
    {
        apply_mass_diag_true(x_true, y_true);
    }

    void apply_minv_half_true(const mfem::Vector &x_true, mfem::Vector &y_true) const;
    void ApplyMinvHalfTrue(const mfem::Vector &x_true, mfem::Vector &y_true) const
    {
        apply_minv_half_true(x_true, y_true);
    }

  private:
    static constexpr int kSpatialDimension = 3;

    static dft::DFTMesh *require_mesh_(const std::shared_ptr<dft::DFTMesh> &mesh);
    static mfem::IntegrationRule make_hex_tensor_gll_rule_(int order);

    void build_diagonal_mass_and_minv_half_();

  private:
    int m_order;
    int m_vdim;

    std::shared_ptr<dft::DFTMesh> m_mesh_owner;
    dft::DFTMesh *m_dft_mesh{nullptr};
    mfem::Mesh *m_mesh{nullptr};
    mfem::H1_FECollection m_fec;
    mfem::FiniteElementSpace m_fes;

    mfem::Vector m_mass_diag_true;
    mfem::Vector m_minv_half_diag_true;
};

#endif // DFT_FESPACE_H
