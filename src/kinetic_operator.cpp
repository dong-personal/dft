#include "kinetic_operator.hpp"

#include <iostream>

KineticOperator::KineticOperator(const DFTGLLHexSpace &space)
    : m_space(space),
      m_fes(const_cast<mfem::FiniteElementSpace *>(&space.Space())),
      m_order(space.Order())
{
}

void KineticOperator::Assemble()
{
    MFEM_VERIFY(m_space.Mesh().GetElementGeometry(0) == mfem::Geometry::CUBE,
                "Expected HEX (Geometry::CUBE)");

    mfem::BilinearForm kinetic_form(m_fes);
    auto *diffusion_integrator = new mfem::DiffusionIntegrator(m_half_coefficient);

    mfem::IntegrationRule integration_rule = make_tensor_gll_rule_(m_order);
    diffusion_integrator->SetIntRule(&integration_rule);

    kinetic_form.AddDomainIntegrator(diffusion_integrator);
    kinetic_form.Assemble();
    kinetic_form.Finalize();

    mfem::Array<int> essential_true_dofs;
    mfem::OperatorHandle matrix_handle;
    kinetic_form.FormSystemMatrix(essential_true_dofs, matrix_handle);

    auto *matrix_ptr = matrix_handle.Is<mfem::SparseMatrix>();
    MFEM_VERIFY(matrix_ptr != nullptr, "Expected FormSystemMatrix to produce a SparseMatrix");

    mfem::SparseMatrix matrix_copy(*matrix_ptr);
    m_matrix.Swap(matrix_copy);

    MFEM_VERIFY(m_matrix.Height() > 0, "Kinetic matrix has zero height");
    MFEM_VERIFY(m_matrix.NumNonZeroElems() > 0, "Kinetic matrix assembly produced zero nnz");

    m_is_assembled = true;
}

bool KineticOperator::IsAssembled() const
{
    return m_is_assembled;
}

const mfem::SparseMatrix &KineticOperator::Matrix() const
{
    MFEM_VERIFY(m_is_assembled, "Kinetic operator has not been assembled");
    return m_matrix;
}

void KineticOperator::Mult(const mfem::Vector &x_true, mfem::Vector &y_true) const
{
    MFEM_VERIFY(m_is_assembled, "Kinetic operator has not been assembled");
    MFEM_VERIFY(x_true.Size() == m_matrix.Width(), "Input size must match kinetic operator width");

    y_true.SetSize(m_matrix.Height());
    m_matrix.Mult(x_true, y_true);
}

double KineticOperator::Energy(const mfem::Vector &x_true) const
{
    mfem::Vector hx;
    Mult(x_true, hx);
    return x_true * hx;
}

void KineticOperator::PrintInfo(std::ostream &os) const
{
    MFEM_VERIFY(m_is_assembled, "Kinetic operator has not been assembled");

    os << "Kinetic matrix size: " << m_matrix.Height() << " x " << m_matrix.Width() << '\n';
    os << "Kinetic matrix nnz: " << m_matrix.NumNonZeroElems() << '\n';
}

mfem::IntegrationRule KineticOperator::make_tensor_gll_rule_(int order)
{
    const int points_per_dimension = order + 1;

    mfem::IntegrationRule line_rule(points_per_dimension);
    mfem::QuadratureFunctions1D quadrature_functions;
    quadrature_functions.GaussLobatto(points_per_dimension, &line_rule);

    return mfem::IntegrationRule(line_rule, line_rule, line_rule);
}
