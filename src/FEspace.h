#pragma once
#include "DFTMesh.h"
#include "mfem.hpp"
#include "mfem/mesh/mesh.hpp"
#include <cmath>
#include <memory>
#include <vector>

class DFTGLLHexSpace
{
  public:
    explicit DFTGLLHexSpace(std::shared_ptr<dft::DFTMesh> mesh, int order, int vdim = 1);
    DFTGLLHexSpace(dft::DFTMesh &mesh, int order, int vdim = 1);

    DFTGLLHexSpace(const DFTGLLHexSpace &) = delete;
    DFTGLLHexSpace &operator=(const DFTGLLHexSpace &) = delete;
    DFTGLLHexSpace(DFTGLLHexSpace &&) = delete;
    DFTGLLHexSpace &operator=(DFTGLLHexSpace &&) = delete;

    mfem::Mesh &Mesh() { return *mesh_; }
    const mfem::Mesh &Mesh() const { return *mesh_; }

    dft::DFTMesh &MeshSource() { return *dft_mesh_; }
    const dft::DFTMesh &MeshSource() const { return *dft_mesh_; }

    mfem::FiniteElementSpace &Space() { return fes_; }
    const mfem::FiniteElementSpace &Space() const { return fes_; }

    int Order() const { return order_; }
    int getDOF() const { return fes_.GetTrueVSize(); } // 真实自由度数（周期性处理后）

    const mfem::Vector &MassDiagTrue() const { return Mdiag_true_; }
    const mfem::Vector &MinvHalfDiagTrue() const { return Minvhalf_true_; }
    std::vector<dft::DFTMesh::Vec3> TrueDofCoordinates() const;

    void ApplyMassDiagTrue(const mfem::Vector &x_true, mfem::Vector &y_true) const;
    void ApplyMinvHalfTrue(const mfem::Vector &x_true, mfem::Vector &y_true) const;

  private:
    static dft::DFTMesh *RequireMesh_(const std::shared_ptr<dft::DFTMesh> &mesh);

    // ---------------- data ----------------
    int order_;
    int vdim_;
    // mfem::Ordering::Type ordering_;

    std::shared_ptr<dft::DFTMesh> mesh_owner_;
    dft::DFTMesh *dft_mesh_{nullptr};
    mfem::Mesh *mesh_{nullptr};
    mfem::H1_FECollection fec_;
    mfem::FiniteElementSpace fes_;

    mfem::Vector Mdiag_true_;
    mfem::Vector Minvhalf_true_;

  private:
    // ---------- tensor GLL rule ----------
    // 定义积分规则, 这里的积分节点与有限元基函数节点相同（GLL节点）, 但是不同的东西, 可以使用更高阶的规则来提高积分精度
    static mfem::IntegrationRule MakeHexTensorGLLRule_(int p)
    {
        const int np = p + 1;

        mfem::IntegrationRule ir1(np);
        mfem::QuadratureFunctions1D qf;
        qf.GaussLobatto(np, &ir1);

        return mfem::IntegrationRule(ir1, ir1, ir1);
    }

    // ---------- diagonal mass + Minvhalf ----------
    // 计算质量矩阵的对角元素和其逆平方根
    void BuildDiagonalMassAndMinvhalf_();
};
