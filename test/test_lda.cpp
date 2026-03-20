#include "mfem.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

using namespace mfem;

// ------------------------ POSCAR 读取 ------------------------
struct Atom {
  int Z;
  Vector r;
};
struct Cell {
  DenseMatrix A; // 3x3 列为晶格矢量 a1 a2 a3（Cartesian）
  std::vector<Atom> atoms;
  int nelec = 0;
};

static Vector ParseVec3(std::istream &in) {
  Vector v(3);
  in >> v[0] >> v[1] >> v[2];
  return v;
}

// 极简 POSCAR：支持 VASP5（元素名 + 个数）与 direct/cart
static Cell ReadPOSCAR(const std::string &path,
                       const std::map<std::string, int> &Zmap,
                       const std::map<std::string, int> &Zval_map) {
  std::ifstream f(path);
  if (!f) {
    throw std::runtime_error("Cannot open POSCAR.");
  }

  Cell cell;
  cell.A.SetSize(3, 3);

  std::string line;
  std::getline(f, line); // comment
  double scale;
  f >> scale;

  Vector a1 = ParseVec3(f), a2 = ParseVec3(f), a3 = ParseVec3(f);
  a1 *= scale;
  a2 *= scale;
  a3 *= scale;

  // A columns
  for (int i = 0; i < 3; i++) {
    cell.A(i, 0) = a1[i];
    cell.A(i, 1) = a2[i];
    cell.A(i, 2) = a3[i];
  }

  std::getline(f, line); // end of line
  std::getline(f, line); // element symbols OR counts

  std::stringstream ss(line);
  std::vector<std::string> symbols;
  std::string tok;
  while (ss >> tok)
    symbols.push_back(tok);

  std::getline(f, line);
  std::stringstream ss2(line);
  std::vector<int> counts;
  int c;
  while (ss2 >> c)
    counts.push_back(c);

  bool has_symbols =
      !symbols.empty() && (int)symbols.size() == (int)counts.size();
  if (!has_symbols) {
    // VASP4 style: line we read as symbols is actually counts
    // shift: symbols unknown -> treat as one species
    counts.clear();
    std::stringstream ss3;
    ss3.str(line);
    while (ss3 >> c)
      counts.push_back(c);
    symbols.assign(counts.size(), "X");
  }

  std::getline(f, line); // possibly "Selective dynamics" or coord type
  bool selective = false;
  if (line.size() && (line[0] == 'S' || line[0] == 's')) {
    selective = true;
    std::getline(f, line);
  }

  bool direct = true;
  if (line.size() &&
      (line[0] == 'C' || line[0] == 'c' || line[0] == 'K' || line[0] == 'k'))
    direct = false;

  // build atoms
  cell.atoms.clear();
  for (size_t s = 0; s < counts.size(); s++) {
    int Z = 1;
    if (symbols[s] != "X") {
      auto it = Zmap.find(symbols[s]);
      if (it != Zmap.end())
        Z = it->second;
    }
    int Zval = Z; // fallback
    auto itv = Zval_map.find(symbols[s]);
    if (itv != Zval_map.end())
      Zval = itv->second;

    for (int n = 0; n < counts[s]; n++) {
      Vector x = ParseVec3(f);
      if (selective) {
        std::string t1, t2, t3;
        f >> t1 >> t2 >> t3;
      }
      Atom a;
      a.Z = Zval; // 用“价电子数”当离子电荷(最简模型)
      a.r.SetSize(3);

      if (direct) {
        // r = A * x (x为分数坐标)
        cell.A.Mult(x, a.r);
      } else {
        a.r = x; // already Cartesian
      }
      cell.atoms.push_back(a);
      cell.nelec += a.Z; // 总电子数（最简）
    }
  }
  return cell;
}

// ------------------------ GLL 谱元空间 + 对角质量 ------------------------
class GLLHexSpace {
public:
  GLLHexSpace(const Cell &cell, int nx, int ny, int nz, int order, int vdim = 1,
              Ordering::Type ordering = Ordering::byNODES, double tol = 1e-12)
      : p_(order), vdim_(vdim), ordering_(ordering),
        mesh_(BuildPeriodicHexMesh_(cell, nx, ny, nz, tol)),
        fec_(order, 3, BasisType::GaussLobatto),
        fes_(&mesh_, &fec_, vdim_, ordering_) {
    BuildDiagonalMassAndMinvhalf_();
  }

  Mesh &mesh() { return mesh_; }
  FiniteElementSpace &fes() { return fes_; }
  const Vector &MinvHalfDiagTrue() const { return Minvhalf_true_; }
  const Vector &MassDiagTrue() const { return Mdiag_true_; }

  void ApplyMinvHalf(const Vector &x_true, Vector &y_true) const {
    y_true.SetSize(x_true.Size());
    for (int i = 0; i < x_true.Size(); i++)
      y_true(i) = Minvhalf_true_(i) * x_true(i);
  }

private:
  int p_, vdim_;
  Ordering::Type ordering_;
  Mesh mesh_;
  H1_FECollection fec_;
  FiniteElementSpace fes_;
  Vector Mdiag_true_, Minvhalf_true_;

  static Mesh BuildPeriodicHexMesh_(const Cell &cell, int nx, int ny, int nz,
                                    double tol) {
    Mesh m0 =
        Mesh::MakeCartesian3D(nx, ny, nz, Element::HEXAHEDRON, 1.0, 1.0, 1.0);

    // 把单位立方体 [0,1]^3 映射到晶胞： r = A * x
    // 注意 m0 的顶点坐标在 [0,1]（MakeCartesian3D）
    for (int vi = 0; vi < m0.GetNV(); vi++) {
      double *v = m0.GetVertex(vi);
      Vector x(3);
      x[0] = v[0];
      x[1] = v[1];
      x[2] = v[2];
      Vector r(3);
      cell.A.Mult(x, r);
      v[0] = r[0];
      v[1] = r[1];
      v[2] = r[2];
    }

    // 周期映射：沿 a1,a2,a3 方向平移一个晶格矢量
    std::vector<Vector> translations;
    translations.reserve(3);
    Vector t1(3), t2(3), t3(3);
    for (int i = 0; i < 3; i++) {
      t1[i] = cell.A(i, 0);
      t2[i] = cell.A(i, 1);
      t3[i] = cell.A(i, 2);
    }
    translations.push_back(t1);
    translations.push_back(t2);
    translations.push_back(t3);

    std::vector<int> vmap = m0.CreatePeriodicVertexMapping(translations, tol);
    Mesh mp = Mesh::MakePeriodic(m0, vmap);
    return mp;
  }

  static IntegrationRule MakeHexTensorGLLRule_(int p) {
    int np = p + 1;
    IntegrationRule ir1(np);
    QuadratureFunctions1D qf;
    qf.GaussLobatto(np, &ir1);
    return IntegrationRule(ir1, ir1, ir1);
  }

  void BuildDiagonalMassAndMinvhalf_() {
    BilinearForm m(&fes_);
    auto *mi = new MassIntegrator();
    IntegrationRule ir = MakeHexTensorGLLRule_(p_);
    mi->SetIntRule(&ir);
    m.AddDomainIntegrator(mi);
    m.Assemble();
    m.Finalize();

    SparseMatrix M;
    Array<int> ess;
    m.FormSystemMatrix(ess, M);

    Mdiag_true_.SetSize(M.Height());
    M.GetDiag(Mdiag_true_);
    Minvhalf_true_.SetSize(Mdiag_true_.Size());
    for (int i = 0; i < Mdiag_true_.Size(); i++)
      Minvhalf_true_(i) = 1.0 / std::sqrt(Mdiag_true_(i));
  }
};

// ------------------------ LDA: 交换(最简 Slater) ------------------------
// n: 电子密度 (a.u.)，返回 Vx 与 eps_x（每电子交换能）
static inline void LDA_X_Slater(double n, double &Vx, double &ex) {
  if (n <= 1e-14) {
    Vx = 0.0;
    ex = 0.0;
    return;
  }
  const double c = std::pow(3.0 / M_PI, 1.0 / 3.0);
  ex = -0.75 * c * std::pow(n, 1.0 / 3.0); // eps_x(n)
  Vx = -c * std::pow(n, 1.0 / 3.0);        // d(n*eps_x)/dn
}

// ------------------------ 离子局域势：高斯核“赝势”模型
// ------------------------ V_ion(r) = sum_I (-Z_I) * erf(|r-RI|/sigma)/|r-RI|
// (避免奇点) 这里用一个简单近似：-Z * exp(-r^2/2σ^2)/sqrt(r^2+ε) 也行
struct IonicPotentialCoeff : public Coefficient {
  const Cell &cell;
  double sigma;
  IonicPotentialCoeff(const Cell &c, double s) : cell(c), sigma(s) {}
  virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip) {
    Vector x;
    T.Transform(ip, x);
    double V = 0.0;
    for (auto &a : cell.atoms) {
      Vector d = x;
      d -= a.r;
      double r2 = d * d;
      double r = std::sqrt(r2 + 1e-12);
      double g = std::exp(-0.5 * r2 / (sigma * sigma));
      V += -a.Z * g / r; // 你后续替换成真实赝势即可
    }
    return V;
  }
};

// ------------------------ 密度：rho = 2 * sum_{i=1}^{Nocc} |psi_i|^2
// (Γ点自旋简并) ------------------------
static void BuildDensity(const FiniteElementSpace &fes,
                         const std::vector<GridFunction *> &psis, int nocc,
                         GridFunction &rho) {
  rho = 0.0;
  for (int i = 0; i < nocc; i++) {
    // rho += 2 * psi^2
    GridFunction &p = *psis[i];
    for (int j = 0; j < rho.Size(); j++)
      rho(j) += 2.0 * p(j) * p(j);
  }
}

// ------------------------ 主程序：SCF 骨架 ------------------------
int main(int argc, char *argv[]) {
  //   Device device("cpu");
  //   device.Print();

  OptionsParser args(argc, argv);
  const char *poscar = "POSCAR";
  int nx = 12, ny = 12, nz = 12;
  int p = 6;
  int nstates = 8;
  int scf_max = 50;
  double mix = 0.3;
  double sigma_ion = 0.15;

  args.AddOption(&poscar, "-p", "--poscar", "POSCAR file.");
  args.AddOption(&nx, "-nx", "--nx", "mesh nx");
  args.AddOption(&ny, "-ny", "--ny", "mesh ny");
  args.AddOption(&nz, "-nz", "--nz", "mesh nz");
  args.AddOption(&p, "-o", "--order", "GLL order p");
  args.AddOption(&nstates, "-n", "--nstates", "number of KS states");
  args.AddOption(&scf_max, "-s", "--scf", "scf max iters");
  args.AddOption(&mix, "-m", "--mix", "density mixing");
  args.AddOption(&sigma_ion, "-si", "--sigma-ion", "ionic smoothing");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cout);
    return 1;
  }
  args.PrintOptions(std::cout);

  // 元素/价电子映射（你后续可接真实赝势/价电子表）
  std::map<std::string, int> Zmap = {{"H", 1}, {"C", 6},   {"N", 7},
                                     {"O", 8}, {"Si", 14}, {"Al", 13}};
  std::map<std::string, int> Zval = {{"H", 1}, {"C", 4},  {"N", 5},
                                     {"O", 6}, {"Si", 4}, {"Al", 3}};

  Cell cell = ReadPOSCAR(poscar, Zmap, Zval);
  std::cout << "Atoms: " << cell.atoms.size() << "  nelec=" << cell.nelec
            << "\n";

  // 1) 空间：GLL HEX 周期
  GLLHexSpace space(cell, nx, ny, nz, p, /*vdim=*/1, Ordering::byNODES);
  auto &fes = space.fes();

  // 2) 初始轨道
  std::vector<GridFunction *> psi(nstates);
  for (int i = 0; i < nstates; i++) {
    psi[i] = new GridFunction(&fes);
    psi[i]->Randomize();
  }

  // 3) ρ, Veff
  GridFunction rho(&fes), rho_new(&fes);
  GridFunction Veff_gf(&fes);

  // 4) 刚度矩阵（动能）与质量对角
  BilinearForm a(&fes);
  a.AddDomainIntegrator(new DiffusionIntegrator); // ∫ ∇φ·∇ψ
  a.Assemble();
  a.Finalize();
  SparseMatrix K;
  Array<int> ess;
  a.FormSystemMatrix(ess, K);

  // 5) 离子势系数
  IonicPotentialCoeff Vion(cell, sigma_ion);

  // 6) Poisson: -Δ VH = 4π ρ  （周期下要处理零均值；这里做最简：强制ρ零均值）
  //    你后续建议用 FFT/Ewald 或加约束消除零模。
  BilinearForm lap(&fes);
  lap.AddDomainIntegrator(new DiffusionIntegrator);
  lap.Assemble();
  lap.Finalize();
  SparseMatrix L;
  lap.FormSystemMatrix(ess, L);
  CGSolver cg;
  cg.SetRelTol(1e-10);
  cg.SetAbsTol(0.0);
  cg.SetMaxIter(500);
  cg.SetPrintLevel(0);
  DSmoother prec(L);
  cg.SetPreconditioner(prec);
  cg.SetOperator(L);

  // 7) SCF
  int nocc = cell.nelec / 2; // 自旋简并
  if (nocc > nstates) {
    std::cerr << "nstates too small for occupations.\n";
    return 2;
  }

  rho = 0.0; // 初始密度
  for (int it = 0; it < scf_max; it++) {
    // (a) 计算 VH：解 -Δ VH = 4π (ρ - <ρ>)
    double mean = rho.Sum() / rho.Size();
    GridFunction rhs(&fes);
    rhs = 0.0;
    for (int i = 0; i < rhs.Size(); i++)
      rhs(i) = 4.0 * M_PI * (rho(i) - mean);

    GridFunction VH(&fes);
    VH = 0.0;
    Vector b, x;
    rhs.GetTrueDofs(b);
    VH.GetTrueDofs(x);
    cg.Mult(b, x);
    VH.SetFromTrueDofs(x);

    // (b) 计算 Vxc（最简 Slater X）与 Vion
    for (int j = 0; j < Veff_gf.Size(); j++) {
      double n = std::max(rho(j), 0.0);
      double Vx, ex;
      LDA_X_Slater(n, Vx, ex);

      // Vion 需要在节点处评估：用一个系数投影（最简：用
      // GridFunctionCoefficient） 这里更简单：用 ProjectCoefficient
      // 每步投影一次（成本高但最小可用）
      Veff_gf(j) = Vx; // 先放 Vxc
    }

    GridFunction Vion_gf(&fes);
    Vion_gf.ProjectCoefficient(Vion);
    for (int j = 0; j < Veff_gf.Size(); j++)
      Veff_gf(j) += VH(j) + Vion_gf(j);

    // (c) 组装 Hamiltonian：H = 1/2 K + Veff * M  (在 GLL 同位下势能近似为对角)
    //     这里采取“矩阵形式”最小实现：H = 0.5*K + diag(Veff)*Mdiag （true
    //     dof上） 简化：直接做标准本征：Htilde = Minvhalf * H * Minvhalf
    //     你后续可以改成 matrix-free + LOBPCG。

    // 建一个 Operator 来做 y = Htilde x
    const Vector &Minv = space.MinvHalfDiagTrue();
    const Vector &Mdiag = space.MassDiagTrue();

    // 先把 Veff 投影到 true dof（节点型空间通常相同）
    Vector Veff_true;
    Veff_gf.GetTrueDofs(Veff_true);

    auto Htilde_mult = [&](const Vector &x_in, Vector &y_out) {
      // tmp = Minv * x
      Vector tmp(x_in.Size());
      for (int i = 0; i < x_in.Size(); i++)
        tmp(i) = Minv(i) * x_in(i);

      // y = 0.5*K*tmp + (Veff* Mdiag) * tmp
      Vector y(tmp.Size());
      y = 0.0;
      K.Mult(tmp, y);
      y *= 0.5;
      for (int i = 0; i < y.Size(); i++)
        y(i) += (Veff_true(i) * Mdiag(i)) * tmp(i);

      // y_out = Minv * y
      y_out.SetSize(y.Size());
      for (int i = 0; i < y.Size(); i++)
        y_out(i) = Minv(i) * y(i);
    };

    // (d) 求最低 nstates 个本征对（最小实现：用 MFEM 的 HypreLOBPCG
    // 需要并行/Par）
    //     这里给“骨架”：你可以接 MFEM 的 eigensolver 或外部 Spectra/ARPACK。
    //     ——为了让代码不炸，这里用“伪迭代”：对每个轨道做几步幂迭代/预条件（仅示意）。
    //     你要真正 DFT，本行必须换成正规本征求解器。

    for (int s = 0; s < nstates; s++) {
      Vector x;
      psi[s]->GetTrueDofs(x);
      for (int kstep = 0; kstep < 10; kstep++) {
        Vector y;
        Htilde_mult(x, y);
        double nrm = y.Norml2();
        if (nrm > 0)
          y /= nrm;
        x = y;
      }
      // 回填
      psi[s]->SetFromTrueDofs(x);
    }

    // (e) 用新轨道构造新密度 + mixing
    BuildDensity(fes, psi, nocc, rho_new);

    // 线性混合
    for (int j = 0; j < rho.Size(); j++)
      rho(j) = (1.0 - mix) * rho(j) + mix * rho_new(j);

    // (f) 收敛判据
    Vector dr(rho.Size());
    for (int j = 0; j < rho.Size(); j++)
      dr(j) = rho_new(j) - rho(j);
    double res = dr.Norml2() / std::sqrt((double)rho.Size());
    std::cout << "SCF " << it << "  ||drho||=" << res << "\n";
    if (res < 1e-6)
      break;
  }

  for (auto p : psi)
    delete p;
  return 0;
}
