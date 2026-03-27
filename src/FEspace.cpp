#include "FEspace.h"

#include <iostream>
#include <utility>

dft::DFTMesh *DFTGLLHexSpace::RequireMesh_(const std::shared_ptr<dft::DFTMesh> &mesh)
{
    MFEM_VERIFY(mesh != nullptr, "DFTGLLHexSpace requires a non-null DFTMesh shared_ptr");
    return mesh.get();
}

DFTGLLHexSpace::DFTGLLHexSpace(std::shared_ptr<dft::DFTMesh> mesh, int order, int vdim)
    : DFTGLLHexSpace(*RequireMesh_(mesh), order, vdim)
{
    mesh_owner_ = std::move(mesh);
}

DFTGLLHexSpace::DFTGLLHexSpace(dft::DFTMesh &mesh, int order, int vdim)
    : order_(order), // order是有限元基函数的阶数
      vdim_(vdim),   // vdim=n表示每个节点有n个分量（矢量场），也可以是1（标量场）
      mesh_owner_(), dft_mesh_(&mesh), mesh_(&mesh.mesh()),
      fec_(order, 3, mfem::BasisType::GaussLobatto), // 3是空间维度
      fes_(mesh_, &fec_, vdim_)                      //
{
    BuildDiagonalMassAndMinvhalf_();
}

std::vector<dft::DFTMesh::Vec3> DFTGLLHexSpace::TrueDofCoordinates() const
{
    mfem::FiniteElementSpace coord_fes(mesh_, &fec_, 3, mfem::Ordering::byNODES);
    mfem::GridFunction xyz(&coord_fes);

    auto identity = [](const mfem::Vector &x, mfem::Vector &y) { y = x; };
    mfem::VectorFunctionCoefficient xyz_coeff(3, identity);
    xyz.ProjectCoefficient(xyz_coeff);

    mfem::Vector true_xyz;
    xyz.GetTrueDofs(true_xyz);

    MFEM_VERIFY(true_xyz.Size() == 3 * getDOF(), "True coordinate vector size must match 3 * scalar true dofs");

    std::vector<dft::DFTMesh::Vec3> coords(static_cast<std::size_t>(getDOF()));
    for (int i = 0; i < getDOF(); ++i)
    {
        coords[static_cast<std::size_t>(i)] = {true_xyz(3 * i + 0), true_xyz(3 * i + 1), true_xyz(3 * i + 2)};
    }

    return coords;
}

void DFTGLLHexSpace::ApplyMassDiagTrue(const mfem::Vector &x_true, mfem::Vector &y_true) const
{
    MFEM_VERIFY(x_true.Size() == Mdiag_true_.Size(), "Input size must match true-dof mass diagonal size");
    y_true.SetSize(x_true.Size());
    for (int i = 0; i < x_true.Size(); ++i)
    {
        y_true(i) = Mdiag_true_(i) * x_true(i);
    }
}

void DFTGLLHexSpace::ApplyMinvHalfTrue(const mfem::Vector &x_true, mfem::Vector &y_true) const
{
    MFEM_VERIFY(x_true.Size() == Minvhalf_true_.Size(),
                "Input size must match true-dof inverse sqrt mass diagonal size");
    y_true.SetSize(x_true.Size());
    for (int i = 0; i < x_true.Size(); ++i)
    {
        y_true(i) = Minvhalf_true_(i) * x_true(i);
    }
}

void DFTGLLHexSpace::BuildDiagonalMassAndMinvhalf_()
{
    MFEM_VERIFY(mesh_->GetElementGeometry(0) == mfem::Geometry::CUBE, "Expected HEX (Geometry::CUBE)");

    mfem::BilinearForm m(&fes_);
    auto *mi = new mfem::MassIntegrator();

    mfem::IntegrationRule ir = MakeHexTensorGLLRule_(order_);
    mi->SetIntRule(&ir);

    m.AddDomainIntegrator(mi);
    m.Assemble();
    m.Finalize();

    mfem::SparseMatrix M;
    mfem::Array<int> ess_tdof;
    m.FormSystemMatrix(ess_tdof, M);

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
                off_diag_norm += val * val;
                off_diag_count++;
            }
            else if (val < 0.0)
            {
                diag_count += 1;
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

    std::cout << "\n===== Mesh Info =====" << std::endl;
    std::cout << "Elements: " << mesh_->GetNE() << std::endl;
    std::cout << "Vertices: " << mesh_->GetNV() << std::endl;
    std::cout << "DOFs (before periodicity): " << fes_.GetVSize() << std::endl;
    std::cout << "DOFs (after periodicity): " << fes_.GetTrueVSize() << std::endl;
#endif

    Mdiag_true_.SetSize(M.Height());
    M.GetDiag(Mdiag_true_);

    Minvhalf_true_.SetSize(Mdiag_true_.Size());
    for (int i = 0; i < Mdiag_true_.Size(); ++i)
    {
        MFEM_VERIFY(Mdiag_true_(i) > 0.0, "Mass diagonal must be positive");
        Minvhalf_true_(i) = 1.0 / std::sqrt(Mdiag_true_(i));
    }
}
