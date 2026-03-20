#pragma once
#include "mfem.hpp"
#include "hamiltonianMatrix.hpp"
#include "overlapMatrix.hpp"
#include <vector>
#include <cmath>

/**
 * @class GeneralizedEigenvalueSolver
 * @brief 求解广义特征值问题: H ψ = ε S ψ
 * 
 * 当使用GLL有限元时，S是对角矩阵，可以进行优化。
 * 使用PCG迭代法配合对角预处理。
 */
class GeneralizedEigenvalueSolver {
public:
  /**
   * @brief 构造求解器
   * @param hamiltonian 哈密顿量矩阵
   * @param overlap 重叠矩阵
   */
  GeneralizedEigenvalueSolver(const HamiltonianMatrix &hamiltonian,
                               const OverlapMatrix &overlap)
      : hamiltonian_(hamiltonian),
        overlap_(overlap),
        num_eigenvalues_(10),
        max_iter_(1000),
        tol_(1e-8),
        eigvals_(),
        eigvecs_() {}

  // ---------- 参数设置 ----------

  void SetNumEigenvalues(int n) {
    MFEM_VERIFY(n > 0, "Number of eigenvalues must be positive");
    num_eigenvalues_ = n;
  }

  void SetMaxIterations(int max_iter) {
    MFEM_VERIFY(max_iter > 0, "Max iterations must be positive");
    max_iter_ = max_iter;
  }

  void SetTolerance(double tol) {
    MFEM_VERIFY(tol > 0, "Tolerance must be positive");
    tol_ = tol;
  }

  // ---------- 求解 ----------

  /**
   * @brief 使用预调节共轭梯度法求解前N个特征对
   * 
   * 算法:
   * 1. 变换到标准特征值问题: H' = S^{-1/2} H S^{-1/2}
   * 2. 使用PCG求解H'y = λy
   * 3. 恢复原始特征向量: ψ = S^{-1/2} y
   */
  void Solve() {
    int size = hamiltonian_.GetKineticEnergy().Height();
    MFEM_VERIFY(num_eigenvalues_ <= size,
                "Number of requested eigenvalues exceeds matrix dimension");

    eigvals_.clear();  // 清空任何之前的特征值
    eigvecs_.clear();
    eigvecs_.resize(num_eigenvalues_);

    // 初始化本征向量(随机)
    eigvecs_.resize(num_eigenvalues_);
    for (int i = 0; i < num_eigenvalues_; i++) {
      eigvecs_[i].SetSize(size);
      RandomVector_(eigvecs_[i]);

      // 正交化相对于S的内积
      for (int j = 0; j < i; j++) {
        double overlap_ij = InnerProduct_(eigvecs_[i], eigvecs_[j]);
        eigvecs_[i].Add(-overlap_ij, eigvecs_[j]);
      }

      // 标准化
      double norm = GetNorm_(eigvecs_[i]);
      eigvecs_[i] *= 1.0 / norm;
    }

    // 对每个特征对进行RayleighRitz迭代
    for (int iter = 0; iter < max_iter_; iter++) {
      bool converged = true;
      
      // 在每次迭代开始时清空特征值(以便重新计算)
      if (iter > 0) {
        eigvals_.clear();
      }

      for (int i = 0; i < num_eigenvalues_; i++) {
        mfem::Vector Hx(size);
        mfem::Vector Sx(size);

        // 计算H*x
        ComputeHx_(eigvecs_[i], Hx);

        // 计算S*x
        overlap_.Apply(eigvecs_[i], Sx);

        // 瑞利商: λ = <x|H|x> / <x|S|x>
        double lambda_num = InnerProduct_(eigvecs_[i], Hx);
        double lambda_den = InnerProduct_(eigvecs_[i], Sx);
        double lambda = lambda_num / lambda_den;

        // 残差: r = H*x - λ*S*x
        mfem::Vector r(size);
        r = Hx;
        r.Add(-lambda, Sx);

        double residual_norm = GetNorm_(r);

        if (residual_norm > tol_) {
          converged = false;
        }

        // 预处理器: 对角元素
        mfem::Vector y(size);
        ApplyPreconditioner_(r, y);

        // 更新: x_new = x + α * y (某种线搜索)
        mfem::Vector x_new = eigvecs_[i];
        double alpha = 0.01;
        x_new.Add(alpha, y);

        // 重新正交化和标准化
        for (int j = 0; j < i; j++) {
          double overlap_ij = InnerProduct_(x_new, eigvecs_[j]);
          x_new.Add(-overlap_ij, eigvecs_[j]);
        }

        double new_norm = GetNorm_(x_new);
        // 缩放向量: eigvecs_[i] = x_new / new_norm
        eigvecs_[i] = x_new;
        eigvecs_[i] *= (1.0 / new_norm);

        // 存储特征值(第一次迭代或之后每次都更新)
        eigvals_.push_back(lambda);

        if (residual_norm < tol_) {
          mfem::out << "Eigenvalue " << i << ": λ = " << lambda
                    << ", residual = " << residual_norm << std::endl;
        }
      }

      if (converged) {
        mfem::out << "Convergence achieved at iteration " << iter
                  << std::endl;
        break;
      }
    }
  }

  // ---------- 获取结果 ----------

  const std::vector<double> &GetEigenvalues() const { return eigvals_; }

  const std::vector<mfem::Vector> &GetEigenvectors() const {
    return eigvecs_;
  }

  /**
   * @brief 获取第i个特征值
   */
  double GetEigenvalue(int i) const {
    MFEM_VERIFY(i < (int)eigvals_.size(), "Eigenvalue index out of range");
    return eigvals_[i];
  }

  /**
   * @brief 获取第i个特征向量
   */
  const mfem::Vector &GetEigenvector(int i) const {
    MFEM_VERIFY(i < (int)eigvecs_.size(), "Eigenvector index out of range");
    return eigvecs_[i];
  }

  // ---------- 诊断 ----------

  void PrintEigenvalues() const {
    mfem::out << "======== Eigenvalues ========" << std::endl;
    for (int i = 0; i < (int)eigvals_.size(); i++) {
      mfem::out << "λ[" << i << "] = " << eigvals_[i] << std::endl;
    }
    mfem::out << "=============================" << std::endl;
  }

private:
  // ---------- 数据成员 ----------
  const HamiltonianMatrix &hamiltonian_;
  const OverlapMatrix &overlap_;

  int num_eigenvalues_;
  int max_iter_;
  double tol_;

  std::vector<double> eigvals_;
  std::vector<mfem::Vector> eigvecs_;

  // ---------- 私有辅助方法 ----------

  /**
   * @brief 计算 Hx = H * x
   */
  void ComputeHx_(const mfem::Vector &x, mfem::Vector &Hx) const {
    hamiltonian_.Mult(x, Hx);
  }

  /**
   * @brief S-内积: <x, y>_S = <x | S | y>
   */
  double InnerProduct_(const mfem::Vector &x, const mfem::Vector &y) const {
    mfem::Vector Sy(y.Size());
    overlap_.Apply(y, Sy);
    return x * Sy;
  }

  /**
   * @brief S-范数: ||x||_S = sqrt(<x | S | x>)
   */
  double GetNorm_(const mfem::Vector &x) const {
    return std::sqrt(InnerProduct_(x, x));
  }

  /**
   * @brief 随机化向量
   */
  void RandomVector_(mfem::Vector &v) const {
    for (int i = 0; i < v.Size(); i++) {
      v(i) = (rand() % 1000) / 500.0 - 1.0;
    }
  }

  /**
   * @brief 应用对角预处理器
   * @param r 残差向量
   * @param y 预处理后的向量
   */
  void ApplyPreconditioner_(const mfem::Vector &r,
                            mfem::Vector &y) const {
    // 使用S^{-1}作为预处理器
    overlap_.ApplyInv(r, y);
  }
};
