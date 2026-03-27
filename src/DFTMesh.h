#pragma once
#ifndef DFT_DFTMESH_H
#define DFT_DFTMESH_H

#include <array>
#include <cmath>
#include <memory>
#include <vector>

#include "Structure.h"
#include "mfem.hpp"

namespace dft
{

class DFTMesh
{
  public:
    using Vec3 = std::array<double, 3>;

    DFTMesh() = default;
    DFTMesh(const DFTMesh &) = delete;
    DFTMesh &operator=(const DFTMesh &) = delete;
    DFTMesh(DFTMesh &&) = delete;
    DFTMesh &operator=(DFTMesh &&) = delete;

    void set_structure(const Structure *s);
    void set_structure(const Structure &s);
    void set_structure(std::shared_ptr<const Structure> s);

    mfem::Mesh &mesh();
    const mfem::Mesh &mesh() const;
    const Structure *structure() const { return structure_; }
    const Structure::Mat3 &lattice() const { return lattice_; }

    // 正确的 save（MFEM Print 接收 ostream）
    void save(const std::string &file) const;

    void saveVTU(const std::string &file) const;

    // ============================================================
    // 生成“三方向周期”的晶胞平行六面体网格（MFEM 4.8）
    // ============================================================
    //
    // lat: Structure::Mat3，约定为行向量 (a1,a2,a3)
    // (u,v,w) in [0,1]^3  --> r = u*a1 + v*a2 + w*a3
    //
    void init_periodic_cell_from_lattice(const Structure::Mat3 &lat, int nx, int ny, int nz,
                                         mfem::Element::Type type = mfem::Element::HEXAHEDRON,
                                         bool sfc_ordering = true);

    void init_periodic_cell_from_structure(int nx, int ny, int nz, mfem::Element::Type type = mfem::Element::HEXAHEDRON,
                                           bool sfc_ordering = true);

    // ============================================================
    // 在原子附近局部细化（保持周期：细化 base_mesh，再重建 periodic）
    // ============================================================
    void refine_near_atoms(double radius, int levels, double shrink = 1.0);

  private:
    void ensure_mesh() const;
    void ensure_base_mesh() const;
    void ensure_structure() const;

    static std::vector<mfem::Vector> make_translations_from_lattice_(const Structure::Mat3 &lat);

    // 根据 base mesh 重建 periodic mesh（核心：保持周期）
    void rebuild_periodic_from_base_();

    // frac -> cart（按你 Mat3 行向量 a1,a2,a3）
    static mfem::Vector frac_to_cart_(const Structure::Mat3 &lat, const Vec3 &f);

    std::vector<mfem::Vector> atom_positions_cartesian_() const;

    static double dist2_(const mfem::Vector &a, const mfem::Vector &b);

    // 计算 x 与原子 a 的最小像距离平方
    double periodic_dist2_(const mfem::Vector &x, const mfem::Vector &a) const;

    // 对 periodic mesh 更稳的中心点：用 ElementTransformation
    // 在参考单元中心求物理坐标
    static mfem::Vector element_center_via_transform_(mfem::Mesh &mesh, int el_id);

  private:
    std::unique_ptr<mfem::Mesh> mesh_;      // 对外：周期 mesh
    std::unique_ptr<mfem::Mesh> base_mesh_; // 内部：非周期 mesh（细化在它上面做）

    std::shared_ptr<const Structure> structure_owner_;
    const Structure *structure_{nullptr};

    Structure::Mat3 lattice_{};              // 保存晶格（可选）
    std::vector<mfem::Vector> translations_; // 保存周期平移向量（重建周期用）
};

} // namespace dft

#endif // DFT_DFTMESH_H
