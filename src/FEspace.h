#pragma once
#include "DFTMesh.h"
#include "mfem.hpp"
#include "mfem/mesh/mesh.hpp"
#include <cmath>

class DFTGLLHexSpace
{
  public:
    DFTGLLHexSpace(dft::DFTMesh &mesh, int order, int vdim = 1)
        : order_(order),                                 // order是有限元基函数的阶数
          vdim_(vdim),                                   // vdim=n表示每个节点有n个分量（矢量场），也可以是1（标量场）
          fec_(order, 3, mfem::BasisType::GaussLobatto), // 3是空间维度
          fes_(&mesh.mesh(), &fec_, vdim_),              //
          mesh_(mesh.mesh())
    {
        BuildDiagonalMassAndMinvhalf_();
    }

    mfem::FiniteElementSpace &Space() { return fes_; }
    const mfem::FiniteElementSpace &Space() const { return fes_; }

    int Order() const { return order_; }
    int getDOF() const { return fes_.GetTrueVSize(); } // 真实自由度数（周期性处理后）

    const mfem::Vector &MassDiagTrue() const { return Mdiag_true_; }

  private:
    // ---------------- data ----------------
    int order_;
    int vdim_;
    // mfem::Ordering::Type ordering_;

    mfem::Mesh mesh_;
    mfem::H1_FECollection fec_;
    mfem::FiniteElementSpace fes_;

    mfem::Vector Mdiag_true_;

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
    void BuildDiagonalMassAndMinvhalf_()
    {
        MFEM_VERIFY(mesh_.GetElementGeometry(0) == mfem::Geometry::CUBE, "Expected HEX (Geometry::CUBE)");

        mfem::BilinearForm m(&fes_);
        auto *mi = new mfem::MassIntegrator();

        mfem::IntegrationRule ir = MakeHexTensorGLLRule_(order_);
        mi->SetIntRule(&ir);

        m.AddDomainIntegrator(mi);
        m.Assemble();
        m.Finalize();

        // mfem::SparseMatrix M;
        // mfem::Array<int> ess_tdof; // periodic Γ-point
        // m.FormSystemMatrix(ess_tdof, M);
        const mfem::SparseMatrix &M = m.SpMat();

#ifdef TEST_COMPILE
        std::cout << "\n===== Overlap Matrix =====" << std::endl;
        std::cout << "Size: " << M.Height() << " x " << M.Width() << std::endl;
        std::cout << "Non-zeros: " << M.NumNonZeroElems() << std::endl;
        double off_diag_norm = 0.0;
        int off_diag_count = 0;
        int diag_count = 0;
        const int *I = M.GetI();
        const int *J = M.GetJ();
        const double *Data = M.GetData();
        for (int i = 0; i < M.Height(); i++)
        {
            for (int k = I[i]; k < I[i + 1]; k++)
            {
                int j = J[k];         // 列号
                double val = Data[k]; // 值
                if (i != j)
                {
                    // std::cout << "M[" << i << "," << j << "] = " << val << std::endl;
                    off_diag_norm += val * val;
                    off_diag_count++;
                }
                else
                {
                    if (val < 0)
                    {
                        // std::cout << "Diagonal M[" << i << "," << j << "] = " << val << std::endl;
                        diag_count += 1;
                    }
                }
            }
        }

        std::cout << "Off-diagonal elements: " << off_diag_count << std::endl;
        std::cout << "Off-diagonal norm: " << std::sqrt(off_diag_norm) << std::endl;
        std::cout << "Diagonal elements: " << diag_count << std::endl;
        std::cout << "Is diagonal: " << (off_diag_norm < 1e-10 ? "YES" : "NO") << std::endl;
        std::cout << "\n===== Integration Rule Info =====" << std::endl;
        std::cout << "Order: " << order_ << std::endl;
        std::cout << "Expected GLL points per dimension: " << (order_ + 1) << std::endl;

        mfem::IntegrationRule ir_debug = MakeHexTensorGLLRule_(order_);
        std::cout << "Actual integration points: " << ir_debug.GetNPoints() << std::endl;

        // 检查网格
        std::cout << "\n===== Mesh Info =====" << std::endl;
        std::cout << "Elements: " << mesh_.GetNE() << std::endl;
        std::cout << "Vertices: " << mesh_.GetNV() << std::endl;
        std::cout << "DOFs (before periodicity): " << fes_.GetVSize() << std::endl;
        std::cout << "DOFs (after periodicity): " << fes_.GetTrueVSize() << std::endl;
#endif

        Mdiag_true_.SetSize(M.Height());
        M.GetDiag(Mdiag_true_);
    }
};
