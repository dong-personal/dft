# 谱有限元DFT求解器实现指南

## 📚 概述

这份文档详细描述了如何在您的DFT有限元框架中实现**使用GLL正交的谱有限元方法**，以得到对角的重叠矩阵。

## 🎯 为什么选择谱有限元+GLL正交？

### 关键优势

1. **对角重叠矩阵** [34]
   - 使用Gauss-Lobatto-Legendre (GLL)点作为基函数节点
   - GLL正交规则 → S_ij = ∫ φ_i φ_j w(x)dx ≈ δ_ij * w_i
   - 消除了传统有限元中的重叠矩阵反演步骤
   - 显著加快了迭代求解器的收敛

2. **谱精度**
   - 代数收敛速度(指数收敛)vs 代数多项式收敛
   - 用更少的自由度获得更高的精度
   - 特别适合于光滑的势场

3. **硬件效率**
   - 张量积结构 → 高度向量化
   - 高计算密度(FLOPs/字节)
   - 在GPU上表现出色

## 📂 项目结构

```
src/
├── fespace.hpp              ✅ GLL有限元空间(已有)
├── overlapMatrix.hpp        ✨ NEW: 对角重叠矩阵
├── hamiltonianMatrix.hpp    ✨ NEW: Kohn-Sham哈密顿量
├── eigenvalueSolver.hpp     ✨ NEW: 特征值求解器
├── spectralDFTSolver.hpp    ✨ NEW: 主求解器类
└── ...

test/
├── test_spectral_dft.cpp    ✨ NEW: 演示程序
└── ...
```

## 🔧 实现步骤

### 第1步: 验证GLL有限元空间

确保 `fespace.hpp` 中的 `DFTGLLHexSpace` 正确:
- ✅ 使用 `mfem::BasisType::GaussLobatto`
- ✅ 使用GLL正交规则 (`QuadratureFunctions1D::GaussLobatto`)
- ✅ 计算对角质量矩阵

```cpp
// fespace.hpp 中已经实现:
DFTGLLHexSpace(int nx, int ny, int nz, int order, ...)
  : fec_(order, 3, mfem::BasisType::GaussLobatto),  // ← GLL基
    fes_(&mesh_, &fec_, ...) { }
```

### 第2步: 组装重叠矩阵

`overlapMatrix.hpp` 处理:
- 使用GLL正交规则组装质量矩阵
- 提取对角元素
- 预计算 S^{-1} 和 S^{-1/2}

**关键代码**:
```cpp
void BuildOverlapMatrix_() {
  mfem::BilinearForm bf(&fes_);
  auto *mi = new mfem::MassIntegrator();
  
  // 使用GLL正交
  mfem::IntegrationRule ir = MakeTensorGLLRule_(p_);
  mi->SetIntRule(&ir);
  
  bf.AddDomainIntegrator(mi);
  bf.Assemble();
  // S_ij ≈ δ_ij * w_i (对角!)
}
```

### 第3步: 组装哈密顿量

`hamiltonianMatrix.hpp` 处理两部分:

**a) 动能项** (刚度矩阵):
```
T_ij = (1/2) ∫ ∇φ_i · ∇φ_j w(x) dx
```

**b) 势能项**:
```
V_ij = ∫ V(r) φ_i(r) φ_j(r) w(x) dx
```

同样使用GLL正交!

### 第4步: 求解特征值问题

`eigenvalueSolver.hpp` 实现广义特征值问题:

```
H ψ = ε S ψ    (S是对角的!)
```

**优化的算法**:
1. 变换到标准问题: H' = S^{-1/2} H S^{-1/2}
2. 求解: H' y = λ y
3. 恢复: ψ = S^{-1/2} y

由于S是对角的:
- S^{-1/2} 的应用 O(n)
- 显著降低计算成本

### 第5步: 自洽场迭代

`spectralDFTSolver.hpp` 协调整体流程:

```cpp
for (scf_iter = 0; scf_iter < max_iter; scf_iter++) {
  // 1. 组装H[ρ]
  hamiltonian_.AssembleKineticEnergy();
  hamiltonian_.AssembleExternalPotential(V_ext);
  hamiltonian_.SetHartreePotential(rho);
  hamiltonian_.AddXCPotential(vxc);
  
  // 2. 求解特征值
  eigensolver_.Solve();
  
  // 3. 更新密度
  UpdateDensity_();
  
  // 4. 检查收敛
  if (DensityError() < tolerance) break;
}
```

## 📊 算法对比

### 传统有限元 (非对角S)
```
对于每次迭代:
- 求解: (H - λS) v = 0
  需要: S^{-1} (O(n³) 或迭代O(n²))
- 总复杂度: O(n³) 或 O(iter × n²)
```

### 谱有限元+GLL (对角S)
```
对于每次迭代:
- 应用: S^{-1/2} v (O(n))
- 总复杂度: O(iter × n) ← 显著加速!
```

## 🔬 GLL性质详解

### 1. GLL节点和权重

对于1D Legendre多项式 P_p(x), GLL点是:
- P_p'(x_i) = 0  (Legendre多项式的导数零点)
- 加上边界点 x_0 = -1, x_p = 1

```
例如 p=2: x = {-1, 0, 1}
         w = {1/3, 4/3, 1/3}
```

### 2. 正交性

GLL正交规则具有精确性:
```
∫_{-1}^{1} P(x) φ_i(x) φ_j(x) dx = δ_ij w_i  (对于多项式P)
```

其中 φ_i 是p阶Legendre多项式的Lagrange基。

### 3. 在3D张量网格中

```
ir_3d = IntegrationRule(ir_1d, ir_1d, ir_1d)
```

自动给出对角的质量矩阵!

## 💡 最佳实践

### 1. 多项式阶数选择
- **p=2-3**: 快速测试 (少量自由度)
- **p=4-6**: 生产计算 (好的精度/成本比)
- **p=8+**: 高精度计算

```cpp
// 示例配置
int nx = 10, ny = 10, nz = 10;  // 网格划分
int p = 4;                       // 谱阶数
SpectralDFTSolver solver(nx, ny, nz, p);
```

自由度数 = (nx × p) × (ny × p) × (nz × p)

### 2. 网格精细化

不要使用太细的网格 (如 nx=100):
- 谱元素具有高分辨率
- p=4 配 nx=10 ≈ p=2 配 nx=20

```cpp
// 好的配置
solver1(10, 10, 10, 4);   // 40K DoFs, 高精度

// 不好的配置
solver2(100, 100, 100, 2); // 8M DoFs, 浪费
```

### 3. 预处理和收敛

使用对角S的优势:
```cpp
// 对角预处理器(非常高效!)
ApplyPreconditioner(r, y):
  y_i = S_diag_inv_(i) * r_i  // O(n)
```

### 4. 周期边界条件

确保周期网格正确创建:
```cpp
// fespace.hpp 中已实现
mfem::Mesh::MakePeriodic(mesh, vertex_mapping);
```

## 🧪 验证步骤

### 测试1: 对角性检查
```cpp
OverlapMatrix overlap(gll_space);
bool is_diag = overlap.IsDiagonal(1e-10);
assert(is_diag);  // 应该通过!
```

### 测试2: 特征值问题
```cpp
// 对于常数势(粒子盒)
// 特征值应该是 π²/(2L²) × (n_x² + n_y² + n_z²)

hamiltonian.AssembleKineticEnergy();
eigensolver.Solve();
auto eigvals = eigensolver.GetEigenvalues();

// 检查第一个特征值
double expected = M_PI * M_PI / 2.0;  // 简化情况
double error = fabs(eigvals[0] - expected);
assert(error < 1e-4);
```

### 测试3: 能量守恒
```cpp
// 在迭代中, 总能量应该单调降低
for (int iter = 0; iter < 10; iter++) {
  solver.SolveSCF();
  double E = solver.ComputeTotalEnergy();
  assert(E <= E_prev + 1e-6);  // 允许数值误差
  E_prev = E;
}
```

## 📖 参考文献

主要参考论文:

1. **谱元素方法基础**
   - Patera, A. T. (1984). "A spectral element method for fluid dynamics: 
     Laminar flow in a channel expansion"

2. **GLL正交和对角质量矩阵** [34]
   - Canuto, C., et al. (2006). "Spectral Methods: Fundamentals in Single 
     Domains", Springer

3. **DFT中的谱方法**
   - Fattebert, J. L., & Gygi, F. (2002). "Density functional theory for 
     efficient ab initio molecular dynamics simulations in the 
     condensed phase"

4. **MFEM库**
   - Kolev, T., & Vassilevski, P. S. (2012). "MFEM: A Modular Finite Element 
     Methods library"

## 🚀 下一步工作

### 短期 (立即实施)
- [ ] 编译并测试overlapMatrix.hpp
- [ ] 验证重叠矩阵确实是对角的
- [ ] 测试简单的特征值问题(粒子盒)
- [ ] 实现CMake集成

### 中期 (1-2周)
- [ ] 添加真实的Hartree势求解器
- [ ] 实现LDA交换相关泛函
- [ ] 完成自洽场循环
- [ ] 添加密度和势的输出

### 长期 (持续改进)
- [ ] GPU加速(CUDA/HIP)
- [ ] 并行化(MPI)
- [ ] 适应性网格精细化
- [ ] 进阶泛函(GGA, 杂交泛函)

## 📝 常见问题

**Q: 为什么一定要用GLL而不是Gauss正交?**

A: Gauss正交给出密集(dense)的质量矩阵。只有在节点与正交点重合时
(即Gauss-Lobatto)才能得到对角矩阵。

**Q: 对角化对精度有影响吗?**

A: 没有!这不是近似,而是精确的性质。S确实是对角的(在机器精度范围内)。

**Q: 如何选择多项式阶数p?**

A: 
- 问题中光滑性越高 → 用更大的p
- 硬件内存/时间有限 → 用较小的p (p=2-3)
- 一般推荐: p=4 是一个很好的平衡点

**Q: 如何处理势场的奇点(如核势)?**

A: 
1. 在原子核处进行局部网格精细化
2. 使用分段定义的势
3. 或使用赝势方法(pseudopotential)避免奇点

## 附录: 完整示例

见 `test/test_spectral_dft.cpp`

```bash
# 编译
cd build
cmake ..
make

# 运行
./test/test_spectral_dft
```

## 版本历史

- v1.0 (2026-03): 初始实现框架
  - 重叠矩阵
  - 哈密顿量组装
  - 特征值求解器
  - 主求解器类

---

**最后更新**: 2026年3月6日  
**联系**: 在项目中提出Issue或Pull Request
