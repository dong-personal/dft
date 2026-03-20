#pragma once
#include "mfem.hpp"
#include "fespace.hpp"
#include "overlapMatrix.hpp"
#include <cmath>

/**
 * @class HamiltonianMatrix
 * @brief DFT Kohn-Sham哈密顿量矩阵
 * 
 * H = -∇²/2 + V_ext + V_hartree + V_xc
 * 
 * 使用GLL有限元和对角重叠矩阵进行优化
 */
class HamiltonianMatrix {
public:
  /**
   * @brief 构造哈密顿量矩阵
   * @param gll_space GLL有限元空间
   */
  explicit HamiltonianMatrix(const DFTGLLHexSpace &gll_space)
      : gll_space_(gll_space),
        fes_(const_cast<mfem::FiniteElementSpace *>(
            &gll_space.Space())),
        p_(gll_space.Order()),
        overlap_(gll_space),
        T_assembled_(false),
        V_assembled_(false),
        coeff_half_(0.5) {
    // T_ 和 V_ 在 Assemble 时由 BilinearForm 创建
  }

  // ---------- 组装各部分 ----------

  /**
   * @brief 组装动能矩阵 T = -∇²/2
   * 
   * 使用谱元素，动能矩阵通过刚度矩阵计算：
   * T_ij = (1/2) ∫ ∇φ_i · ∇φ_j dr
   */
  void AssembleKineticEnergy() {
    MFEM_VERIFY(
        gll_space_.Mesh().GetElementGeometry(0) == mfem::Geometry::CUBE,
        "Expected hexahedral mesh");

    mfem::BilinearForm bf(fes_);
    // DiffusionIntegrator 需要 Coefficient 参数
    auto *di = new mfem::DiffusionIntegrator(coeff_half_);

    // 使用GLL正交规则
    mfem::IntegrationRule ir = MakeTensorGLLRule_(p_);
    di->SetIntRule(&ir);

    bf.AddDomainIntegrator(di);
    bf.Assemble();
    bf.Finalize();

    mfem::Array<int> ess_tdof;
    bf.FormSystemMatrix(ess_tdof, T_);

    T_assembled_ = true;
  }

  /**
   * @brief 组装外部势矩阵 V_ext
   * @param coeff 外部势系数函数
   */
  void AssembleExternalPotential(mfem::Coefficient &coeff) {
    mfem::BilinearForm bf(fes_);
    auto *mi = new mfem::MassIntegrator(coeff);

    mfem::IntegrationRule ir = MakeTensorGLLRule_(p_);
    mi->SetIntRule(&ir);

    bf.AddDomainIntegrator(mi);
    bf.Assemble();
    bf.Finalize();

    mfem::Array<int> ess_tdof;
    bf.FormSystemMatrix(ess_tdof, V_);

    V_assembled_ = true;
  }

  /**
   * @brief 设置Hartree势(离子-电子库伦相互作用)
   * @param rho 电荷密度向量(true-dof)
   * 
   * V_Hartree 通过求解 ∇²φ = -4π ρ
   */
  void SetHartreePotential(const mfem::Vector &rho) {
    // TODO: 实现Poisson求解器
    hartree_potential_ = rho;
  }

  /**
   * @brief 添加交换相关势(LDA, GGA等)
   * @param vxc 交换相关势向量
   * 
   * TODO: 实现通过GridFunction或FunctionCoefficient
   */
  void AddXCPotential(const mfem::Vector &vxc) {
    // 将向量转换为系数并应用
    // mfem::Coefficient vxc_coeff;  // TODO: 创建GridFunction系数
    // AssembleExternalPotential(vxc_coeff);
  }

  // ---------- 获取矩阵 ----------

  const mfem::SparseMatrix &GetKineticEnergy() const {
    MFEM_VERIFY(T_assembled_, "Kinetic energy not assembled");
    return T_;
  }

  const mfem::SparseMatrix &GetPotential() const {
    MFEM_VERIFY(V_assembled_, "Potential not assembled");
    return V_;
  }

  const OverlapMatrix &GetOverlapMatrix() const { return overlap_; }

  /**
   * @brief 获取完整哈密顿量 H = T + V
   * @return 稀疏矩阵(需要复制)
   */
  mfem::SparseMatrix GetTotal() const {
    MFEM_VERIFY(T_assembled_ && V_assembled_,
                "Not all matrix components assembled");
    mfem::SparseMatrix H = T_;
    H += V_;
    return H;
  }

  // ---------- 矩阵-向量乘积 ----------

  /**
   * @brief 计算 y = H * x
   * @param x_true true-dof向量
   * @param y_true 输出向量
   */
  void Mult(const mfem::Vector &x_true, mfem::Vector &y_true) const {
    MFEM_VERIFY(T_assembled_, "Kinetic energy not assembled");
    
    mfem::Vector Tx(T_.Height());
    T_.Mult(x_true, Tx);
    
    if (V_assembled_) {
      mfem::Vector Vx(V_.Height());
      V_.Mult(x_true, Vx);
      Tx += Vx;
    }
    
    y_true = Tx;
  }

  // ---------- 诊断 ----------

  void PrintInfo() const {
    if (T_assembled_) {
      mfem::out << "Kinetic energy matrix: " << T_.Height() << " x "
                << T_.Width() << ", nnz = " << T_.NumNonZeroElems()
                << std::endl;
    }
    if (V_assembled_) {
      mfem::out << "Potential matrix: " << V_.Height() << " x " << V_.Width()
                << ", nnz = " << V_.NumNonZeroElems() << std::endl;
    }
    overlap_.PrintDiagInfo();
  }

private:
  // ---------- 数据成员 ----------
  const DFTGLLHexSpace &gll_space_;
  mfem::FiniteElementSpace *fes_;
  int p_;

  OverlapMatrix overlap_;
  mfem::SparseMatrix T_;  // 动能矩阵
  mfem::SparseMatrix V_;  // 势能矩阵

  mfem::Vector hartree_potential_;

  bool T_assembled_;
  bool V_assembled_;

  // 常数系数 (1/2 用于动能)
  mfem::ConstantCoefficient coeff_half_;

  // ---------- 私有方法 ----------

  /**
   * @brief 创建张量GLL正交规则
   */
  static mfem::IntegrationRule MakeTensorGLLRule_(int p) {
    const int np = p + 1;
    mfem::IntegrationRule ir1(np);
    mfem::QuadratureFunctions1D qf;
    qf.GaussLobatto(np, &ir1);
    return mfem::IntegrationRule(ir1, ir1, ir1);
  }
};
