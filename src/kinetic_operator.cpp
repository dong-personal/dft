#include "kinetic_operator.hpp"

#include <iostream>

KineticOperator::KineticOperator(const DFTGLLHexSpace &space)
    : space_(space), fes_(const_cast<mfem::FiniteElementSpace *>(&space.Space())), order_(space.Order())
{
}

void KineticOperator::Assemble()
{
    MFEM_VERIFY(space_.Mesh().GetElementGeometry(0) == mfem::Geometry::CUBE, "Expected HEX (Geometry::CUBE)");

    mfem::BilinearForm kinetic_form(fes_);
    auto *diffusion = new mfem::DiffusionIntegrator(half_coeff_);

    mfem::IntegrationRule ir = MakeTensorGLLRule_(order_);
    diffusion->SetIntRule(&ir);

    kinetic_form.AddDomainIntegrator(diffusion);
    kinetic_form.Assemble();
    kinetic_form.Finalize();

    mfem::Array<int> ess_tdof;
    mfem::OperatorHandle matrix_handle;
    kinetic_form.FormSystemMatrix(ess_tdof, matrix_handle);

    auto *matrix_ptr = matrix_handle.Is<mfem::SparseMatrix>();
    MFEM_VERIFY(matrix_ptr != nullptr, "Expected FormSystemMatrix to produce a SparseMatrix");

    mfem::SparseMatrix matrix_copy(*matrix_ptr);
    matrix_.Swap(matrix_copy);
    MFEM_VERIFY(matrix_.Height() > 0, "Kinetic matrix has zero height");
    MFEM_VERIFY(matrix_.NumNonZeroElems() > 0, "Kinetic matrix assembly produced zero nnz");
    assembled_ = true;
}

const mfem::SparseMatrix &KineticOperator::Matrix() const
{
    MFEM_VERIFY(assembled_, "Kinetic operator has not been assembled");
    return matrix_;
}

void KineticOperator::Mult(const mfem::Vector &x_true, mfem::Vector &y_true) const
{
    MFEM_VERIFY(assembled_, "Kinetic operator has not been assembled");
    MFEM_VERIFY(x_true.Size() == matrix_.Width(), "Input size must match kinetic operator width");
    y_true.SetSize(matrix_.Height());
    matrix_.Mult(x_true, y_true);
}

double KineticOperator::Energy(const mfem::Vector &x_true) const
{
    mfem::Vector hx;
    Mult(x_true, hx);
    return x_true * hx;
}

void KineticOperator::PrintInfo(std::ostream &os) const
{
    MFEM_VERIFY(assembled_, "Kinetic operator has not been assembled");
    os << "Kinetic matrix size: " << matrix_.Height() << " x " << matrix_.Width() << '\n';
    os << "Kinetic matrix nnz: " << matrix_.NumNonZeroElems() << '\n';
}

mfem::IntegrationRule KineticOperator::MakeTensorGLLRule_(int order)
{
    const int np = order + 1;
    mfem::IntegrationRule ir1(np);
    mfem::QuadratureFunctions1D qf;
    qf.GaussLobatto(np, &ir1);
    return mfem::IntegrationRule(ir1, ir1, ir1);
}
