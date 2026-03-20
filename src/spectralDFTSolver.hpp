#pragma once
#include "mfem.hpp"
#include "fespace.hpp"
#include "overlapMatrix.hpp"
#include "hamiltonianMatrix.hpp"
#include "eigenvalueSolver.hpp"
#include <vector>
#include <cmath>
#include <iostream>

/**
 * @class SpectralDFTSolver
 * @brief 基于谱有限元和GLL正交的DFT求解器
 * 
 * 该类整合以下组件:
 * - GLL有限元空间(谱元素)
 * - 对角重叠矩阵
 * - Kohn-Sham哈密顿量
 * - 特征值求解器
 * - 自洽场迭代
 */
class SpectralDFTSolver {
public:
  /**
   * @brief 构造DFT求解器
   * @param nx, ny, nz 各方向元素数
   * @param order 谱元素多项式阶数
   */
  SpectralDFTSolver(int nx, int ny, int nz, int order)
      : gll_space_(nx, ny, nz, order),
        hamiltonian_(gll_space_),
        overlap_(gll_space_),
        eigensolver_(hamiltonian_, overlap_),
        num_electrons_(0),
        scf_tol_(1e-6),
        scf_max_iter_(100),
        converged_(false) {
    int size = gll_space_.TrueVSize();
    rho_.SetSize(size);
    rho_old_.SetSize(size);
    rho_ = 0.0;
    rho_old_ = 0.0;
  }

  // ---------- 参数设置 ----------

  /**
   * @brief 设置电子数
   */
  void SetNumElectrons(int n) {
    MFEM_VERIFY(n > 0, "Number of electrons must be positive");
    num_electrons_ = n;
    eigensolver_.SetNumEigenvalues(n);
  }

  /**
   * @brief 设置自洽场迭代的收敛容限
   */
  void SetSCFTolerance(double tol) {
    MFEM_VERIFY(tol > 0, "SCF tolerance must be positive");
    scf_tol_ = tol;
  }

  /**
   * @brief 设置最大SCF迭代次数
   */
  void SetSCFMaxIterations(int max_iter) {
    MFEM_VERIFY(max_iter > 0, "Max SCF iterations must be positive");
    scf_max_iter_ = max_iter;
  }

  // ---------- 初始化和求解 ----------

  /**
   * @brief 初始化求解器(设置初始密度)
   */
  void Initialize() {
    // 初始密度为均匀分布
    double volume = ComputeVolume_();
    double initial_density = num_electrons_ / volume;

    for (int i = 0; i < rho_.Size(); i++) {
      rho_(i) = initial_density;
    }
    rho_old_ = rho_;
  }

  /**
   * @brief 执行自洽场(SCF)迭代
   * 
   * 算法:
   * 1. 给定电荷密度ρ
   * 2. 计算有效势 V_eff = V_ext + V_hartree + V_xc
   * 3. 求解Kohn-Sham特征值问题
   * 4. 更新电荷密度
   * 5. 检查收敛性，若未收敛返回步骤2
   */
  void SolveSCF() {
    MFEM_VERIFY(num_electrons_ > 0,
                "Number of electrons not set. Call SetNumElectrons().");

    mfem::out << "========================================" << std::endl;
    mfem::out << "   Spectral Element DFT Solver (GLL)   " << std::endl;
    mfem::out << "========================================" << std::endl;
    mfem::out << "Polynomial order: " << gll_space_.Order() << std::endl;
    mfem::out << "Number of DoFs: " << gll_space_.TrueVSize() << std::endl;
    mfem::out << "Number of electrons: " << num_electrons_ << std::endl;
    mfem::out << "========================================" << std::endl;

    for (int scf_iter = 0; scf_iter < scf_max_iter_; scf_iter++) {
      mfem::out << "\nSCF Iteration " << scf_iter << std::endl;
      mfem::out << "--------------------------------------" << std::endl;

      // 步骤1: 组装哈密顿量
      hamiltonian_.AssembleKineticEnergy();
      // TODO: 根据当前密度添加势能项

      // 步骤2: 求解特征值问题
      eigensolver_.Solve();
      eigensolver_.PrintEigenvalues();

      // 步骤3: 计算新密度
      UpdateDensity_();

      // 步骤4: 检查收敛性
      double density_error = ComputeDensityError_();
      mfem::out << "Density error: " << density_error << std::endl;

      if (density_error < scf_tol_) {
        mfem::out << "\nSCF Convergence achieved at iteration " << scf_iter
                  << std::endl;
        converged_ = true;
        break;
      }

      rho_old_ = rho_;
    }

    if (!converged_) {
      mfem::out << "Warning: SCF did not converge within " << scf_max_iter_
                << " iterations" << std::endl;
    }
  }

  // ---------- 结果获取 ----------

  bool IsConverged() const { return converged_; }

  const mfem::Vector &GetDensity() const { return rho_; }

  const std::vector<double> &GetKohnShamEigenvalues() const {
    return eigensolver_.GetEigenvalues();
  }

  const std::vector<mfem::Vector> &GetKohnShamEigenvectors() const {
    return eigensolver_.GetEigenvectors();
  }

  /**
   * @brief 计算总能量
   * 
   * E_total = T[ρ] + E_ext[ρ] + E_hartree[ρ] + E_xc[ρ] - μ*N
   * 
   * 其中 T[ρ] 是动能，E_ext 是外部势能，E_hartree 是Hartree能，
   *      E_xc 是交换相关能，μ 是化学势，N 是电子数
   */
  double ComputeTotalEnergy() const {
    // TODO: 实现能量计算
    return 0.0;
  }

  /**
   * @brief 计算带隙(HOMO-LUMO间隙)
   */
  double ComputeBandgap() const {
    const auto &eigvals = eigensolver_.GetEigenvalues();
    if (eigvals.size() < num_electrons_ + 1) {
      return 0.0;  // 数据不足
    }
    // HOMO是第num_electrons-1个特征值(0-indexed)
    // LUMO是第num_electrons个特征值
    return eigvals[num_electrons_] - eigvals[num_electrons_ - 1];
  }

  // ---------- 网格和空间访问 ----------

  const DFTGLLHexSpace &GetFESpace() const { return gll_space_; }

  mfem::FiniteElementSpace &GetFESpace() {
    return const_cast<mfem::FiniteElementSpace &>(gll_space_.Space());
  }

private:
  // ---------- 数据成员 ----------
  DFTGLLHexSpace gll_space_;
  HamiltonianMatrix hamiltonian_;
  OverlapMatrix overlap_;
  GeneralizedEigenvalueSolver eigensolver_;

  int num_electrons_;
  double scf_tol_;
  int scf_max_iter_;
  bool converged_;

  mfem::Vector rho_;      // 当前电荷密度
  mfem::Vector rho_old_;  // 前一步的电荷密度

  // ---------- 私有辅助方法 ----------

  /**
   * @brief 计算模拟盒子的体积
   */
  double ComputeVolume_() const {
    // 对于周期晶胞，体积由晶格参数确定
    // 这里假设单位立方体[0,1]^3的体积为1
    return 1.0;
  }

  /**
   * @brief 从占据的轨道更新电荷密度
   * 
   * ρ(r) = 2 * Σ_{i=1}^{N/2} |ψ_i(r)|²
   * 
   * (因子2是自旋的贡献)
   */
  void UpdateDensity_() {
    int size = rho_.Size();
    rho_ = 0.0;

    const auto &eigvecs = eigensolver_.GetEigenvectors();
    int num_occ = (num_electrons_ + 1) / 2;  // 考虑自旋

    // TODO: 将特征向量转换为实空间值
    // 这需要完整的GridFunction和评估基础设施

    mfem::out << "Updated density (simplified version)" << std::endl;
  }

  /**
   * @brief 计算密度的L2误差估计
   */
  double ComputeDensityError_() const {
    mfem::Vector diff = rho_;
    diff -= rho_old_;
    return diff.Norml2();
  }
};
