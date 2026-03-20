#pragma once
#ifndef DFT_DFTMESH_H
#define DFT_DFTMESH_H

#include <array>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <vector>

#include "Structure.h"
#include "mfem.hpp"

namespace dft
{

class DFTMesh
{
  public:
    using Vec3 = std::array<double, 3>;

    void set_structure(const Structure *s) { structure_ = s; }

    mfem::Mesh &mesh()
    {
        ensure_mesh();
        return *mesh_;
    }
    const mfem::Mesh &mesh() const
    {
        ensure_mesh();
        return *mesh_;
    }

    // 正确的 save（MFEM Print 接收 ostream）
    void save(const std::string &file) const
    {
        ensure_mesh();
        std::ofstream ofs(file);
        if (!ofs)
            throw std::runtime_error("DFTMesh::save: cannot open " + file);
        mesh_->Print(ofs);
    }

    void saveVTU(const std::string &file) const
    {
        ensure_mesh();

        std::ofstream ofs(file);
        if (!ofs)
            throw std::runtime_error("DFTMesh::saveGmsh: cannot open " + file);

        ofs << "<?xml version=\"1.0\"?>\n";
        ofs << "<VTKFile type=\"UnstructuredGrid\" "
               "version=\"0.1\" byte_order=\"LittleEndian\">\n";
        ofs << "<UnstructuredGrid>\n";
        mesh_->PrintVTU(ofs);
        ofs << "</Piece>\n";
        ofs << "</UnstructuredGrid>\n";
        ofs << "</VTKFile>\n";
    }

    // ============================================================
    // 生成“三方向周期”的晶胞平行六面体网格（MFEM 4.8）
    // ============================================================
    //
    // lat: Structure::Mat3，约定为行向量 (a1,a2,a3)
    // (u,v,w) in [0,1]^3  --> r = u*a1 + v*a2 + w*a3
    //
    void init_periodic_cell_from_lattice(const Structure::Mat3 &lat, int nx, int ny, int nz,
                                         mfem::Element::Type type = mfem::Element::HEXAHEDRON, bool sfc_ordering = true)
    {
        if (nx < 3 || ny < 3 || nz < 3)
        {
            throw std::runtime_error("Periodic mesh needs >=2 interior vertices per direction. "
                                     "Use nx,ny,nz >= 3 (elements).");
        }

        // 1) 先建单位立方体 [0,1]^3 的结构化网格（非周期）
        mfem::Mesh m = mfem::Mesh::MakeCartesian3D(nx, ny, nz, type, 1.0, 1.0, 1.0, sfc_ordering);

        // 2) 顶点仿射映射到晶胞平行六面体（必须在 MakePeriodic 之前做）
        for (int i = 0; i < m.GetNV(); ++i)
        {
            double *v = m.GetVertex(i); // v=(u,v,w)
            const double u = v[0], vv = v[1], w = v[2];

            const double x = u * lat[0][0] + vv * lat[1][0] + w * lat[2][0];
            const double y = u * lat[0][1] + vv * lat[1][1] + w * lat[2][1];
            const double z = u * lat[0][2] + vv * lat[1][2] + w * lat[2][2];

            v[0] = x;
            v[1] = y;
            v[2] = z;
        }

        // 保存 lattice 与 translations，后续“细化后重建周期”要用
        lattice_ = lat;
        translations_ = make_translations_from_lattice_(lat);

        // 保存非周期 base mesh（细化永远在它上面做）
        base_mesh_ = std::make_unique<mfem::Mesh>(std::move(m));

        // 3) 构造 periodic mesh
        rebuild_periodic_from_base_();
    }

    void init_periodic_cell_from_structure(int nx, int ny, int nz, mfem::Element::Type type = mfem::Element::HEXAHEDRON,
                                           bool sfc_ordering = true)
    {
        ensure_structure();
        init_periodic_cell_from_lattice(structure_->lattice(), nx, ny, nz, type, sfc_ordering);
    }

    // ============================================================
    // 在原子附近局部细化（保持周期：细化 base_mesh，再重建 periodic）
    // ============================================================
    void refine_near_atoms(double radius, int levels, double shrink = 1.0)
    {
        ensure_base_mesh();
        ensure_structure();
        if (radius <= 0.0 || levels <= 0)
            return;

        base_mesh_->EnsureNCMesh();

        // 原子坐标转为笛卡尔（与 base_mesh 一致）
        std::vector<mfem::Vector> atom_pos = atom_positions_cartesian_();

        for (int lvl = 0; lvl < levels; ++lvl)
        {
            const double r2 = radius * radius;

            mfem::Array<int> refine_list;
            refine_list.Reserve(base_mesh_->GetNE() / 8 + 1);

            for (int el = 0; el < base_mesh_->GetNE(); ++el)
            {
                mfem::Vector c = element_center_via_transform_(*base_mesh_, el);

                bool hit = false;
                for (const auto &rp : atom_pos)
                {
                    if (periodic_dist2_(c, rp) < r2)
                    {
                        hit = true;
                        break;
                    }
                }
                if (hit)
                    refine_list.Append(el);
            }

            if (refine_list.Size() == 0)
                break;

            // 对 hex 的局部细化，需要 NCMesh；这里用 GeneralRefinement 即可
            base_mesh_->GeneralRefinement(refine_list);

            radius *= shrink;
            if (radius <= 0.0)
                break;
        }

        // 关键：细化完成后，重建 periodic mesh，周期性不会丢
        rebuild_periodic_from_base_();
    }

  private:
    void ensure_mesh() const
    {
        if (!mesh_)
            throw std::runtime_error("DFTMesh: mesh not initialized");
    }
    void ensure_base_mesh() const
    {
        if (!base_mesh_)
            throw std::runtime_error("DFTMesh: base mesh not initialized");
    }
    void ensure_structure() const
    {
        if (!structure_)
            throw std::runtime_error("DFTMesh: structure not set");
    }

    static std::vector<mfem::Vector> make_translations_from_lattice_(const Structure::Mat3 &lat)
    {
        std::vector<mfem::Vector> translations;
        translations.reserve(3);

        for (int i = 0; i < 3; ++i)
        {
            mfem::Vector t(3);
            t[0] = lat[i][0];
            t[1] = lat[i][1];
            t[2] = lat[i][2];
            translations.push_back(t);
        }
        return translations;
    }

    // 根据 base mesh 重建 periodic mesh（核心：保持周期）
    void rebuild_periodic_from_base_()
    {
        ensure_base_mesh();
        if (translations_.empty())
        {
            throw std::runtime_error("DFTMesh: translations not set (init first)");
        }

        std::vector<int> v2v = base_mesh_->CreatePeriodicVertexMapping(translations_);
        mfem::Mesh pm = mfem::Mesh::MakePeriodic(*base_mesh_, v2v);
        mesh_ = std::make_unique<mfem::Mesh>(std::move(pm));
    }

    // frac -> cart（按你 Mat3 行向量 a1,a2,a3）
    static mfem::Vector frac_to_cart_(const Structure::Mat3 &lat, const Vec3 &f)
    {
        mfem::Vector r(3);
        r[0] = f[0] * lat[0][0] + f[1] * lat[1][0] + f[2] * lat[2][0];
        r[1] = f[0] * lat[0][1] + f[1] * lat[1][1] + f[2] * lat[2][1];
        r[2] = f[0] * lat[0][2] + f[1] * lat[1][2] + f[2] * lat[2][2];
        return r;
    }

    std::vector<mfem::Vector> atom_positions_cartesian_() const
    {
        std::vector<mfem::Vector> pos;
        pos.reserve(structure_->num_atoms());
        const auto &lat = structure_->lattice();

        for (const auto &a : structure_->atoms())
        {
            const auto &p = a.position();
            if (structure_->coord_type() == Structure::CoordType::Fractional)
            {
                pos.push_back(frac_to_cart_(lat, p));
            }
            else
            {
                mfem::Vector r(3);
                r[0] = p[0];
                r[1] = p[1];
                r[2] = p[2];
                pos.push_back(r);
            }
        }
        return pos;
    }

    static double dist2_(const mfem::Vector &a, const mfem::Vector &b)
    {
        const double dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
        return dx * dx + dy * dy + dz * dz;
    }

    // 计算 x 与原子 a 的最小像距离平方
    double periodic_dist2_(const mfem::Vector &x, const mfem::Vector &a) const
    {
        const auto &lat = structure_->lattice();

        // 晶格向量 a1,a2,a3
        mfem::Vector a1(3), a2(3), a3(3);
        for (int d = 0; d < 3; ++d)
        {
            a1[d] = lat[0][d];
            a2[d] = lat[1][d];
            a3[d] = lat[2][d];
        }

        double min_d2 = std::numeric_limits<double>::max();

        // 27 个周期像
        for (int i = -1; i <= 1; ++i)
        {
            for (int j = -1; j <= 1; ++j)
            {
                for (int k = -1; k <= 1; ++k)
                {
                    const double dx = x[0] - (a[0] + i * a1[0] + j * a2[0] + k * a3[0]);
                    const double dy = x[1] - (a[1] + i * a1[1] + j * a2[1] + k * a3[1]);
                    const double dz = x[2] - (a[2] + i * a1[2] + j * a2[2] + k * a3[2]);

                    const double d2 = dx * dx + dy * dy + dz * dz;
                    if (d2 < min_d2)
                        min_d2 = d2;
                }
            }
        }
        return min_d2;
    }

    // 对 periodic mesh 更稳的中心点：用 ElementTransformation
    // 在参考单元中心求物理坐标
    static mfem::Vector element_center_via_transform_(mfem::Mesh &mesh, int el_id)
    {
        mfem::ElementTransformation *T = mesh.GetElementTransformation(el_id);
        mfem::IntegrationPoint ip;
        // 参考单元中心：对常见单元(HEX/TET)用 (0,0,0) 很通用
        ip.Set3(0.5, 0.5, 0.5);

        mfem::Vector x(3);
        T->Transform(ip, x);
        return x;
    }

  private:
    std::unique_ptr<mfem::Mesh> mesh_;      // 对外：周期 mesh
    std::unique_ptr<mfem::Mesh> base_mesh_; // 内部：非周期 mesh（细化在它上面做）

    const Structure *structure_{nullptr};

    Structure::Mat3 lattice_{};              // 保存晶格（可选）
    std::vector<mfem::Vector> translations_; // 保存周期平移向量（重建周期用）
};

} // namespace dft

#endif // DFT_DFTMESH_H
