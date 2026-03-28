#include "fespace.h"

#include <cmath>
#include <iostream>
#include <utility>

dft::DFTMesh *DFTGLLHexSpace::require_mesh_(const std::shared_ptr<dft::DFTMesh> &mesh)
{
    MFEM_VERIFY(mesh != nullptr, "DFTGLLHexSpace requires a non-null DFTMesh shared_ptr");
    return mesh.get();
}

DFTGLLHexSpace::DFTGLLHexSpace(std::shared_ptr<dft::DFTMesh> mesh, int order, int vdim)
    : DFTGLLHexSpace(*require_mesh_(mesh), order, vdim)
{
    m_mesh_owner = std::move(mesh);
}

DFTGLLHexSpace::DFTGLLHexSpace(dft::DFTMesh &mesh, int order, int vdim)
    : m_order(order),
      m_vdim(vdim),
      m_mesh_owner(),
      m_dft_mesh(&mesh),
      m_mesh(&mesh.mesh()),
      m_fec(order, kSpatialDimension, mfem::BasisType::GaussLobatto),
      m_fes(m_mesh, &m_fec, m_vdim)
{
    build_diagonal_mass_and_minv_half_();
}

std::vector<DFTGLLHexSpace::SpatialPoint> DFTGLLHexSpace::true_dof_coordinates() const
{
    mfem::FiniteElementSpace coordinate_space(m_mesh, &m_fec, kSpatialDimension, mfem::Ordering::byNODES);
    mfem::GridFunction coordinates(&coordinate_space);

    auto identity = [](const mfem::Vector &input_point, mfem::Vector &output_point) { output_point = input_point; };
    mfem::VectorFunctionCoefficient coordinate_coefficient(kSpatialDimension, identity);
    coordinates.ProjectCoefficient(coordinate_coefficient);

    mfem::Vector true_coordinates;
    coordinates.GetTrueDofs(true_coordinates);

    MFEM_VERIFY(true_coordinates.Size() == kSpatialDimension * true_dof_count(),
                "True coordinate vector size must match kSpatialDimension * true dof count");

    std::vector<SpatialPoint> points;
    points.reserve(static_cast<std::size_t>(true_dof_count()));
    for (int dof_index = 0; dof_index < true_dof_count(); ++dof_index)
    {
        SpatialPoint point(typename SpatialPoint::ShapeType{3});
        for (int dim = 0; dim < kSpatialDimension; ++dim)
        {
            point[dim] = true_coordinates(kSpatialDimension * dof_index + dim);
        }
        points.push_back(std::move(point));
    }

    return points;
}

void DFTGLLHexSpace::apply_mass_diag_true(const mfem::Vector &x_true, mfem::Vector &y_true) const
{
    MFEM_VERIFY(x_true.Size() == m_mass_diag_true.Size(),
                "Input size must match true-dof mass diagonal size");

    y_true.SetSize(x_true.Size());
    for (int dof_index = 0; dof_index < x_true.Size(); ++dof_index)
    {
        y_true(dof_index) = m_mass_diag_true(dof_index) * x_true(dof_index);
    }
}

void DFTGLLHexSpace::apply_minv_half_true(const mfem::Vector &x_true, mfem::Vector &y_true) const
{
    MFEM_VERIFY(x_true.Size() == m_minv_half_diag_true.Size(),
                "Input size must match true-dof inverse sqrt mass diagonal size");

    y_true.SetSize(x_true.Size());
    for (int dof_index = 0; dof_index < x_true.Size(); ++dof_index)
    {
        y_true(dof_index) = m_minv_half_diag_true(dof_index) * x_true(dof_index);
    }
}

mfem::IntegrationRule DFTGLLHexSpace::make_hex_tensor_gll_rule_(int order)
{
    const int points_per_dimension = order + 1;

    mfem::IntegrationRule line_rule(points_per_dimension);
    mfem::QuadratureFunctions1D quadrature_functions;
    quadrature_functions.GaussLobatto(points_per_dimension, &line_rule);

    return mfem::IntegrationRule(line_rule, line_rule, line_rule);
}

void DFTGLLHexSpace::build_diagonal_mass_and_minv_half_()
{
    MFEM_VERIFY(m_mesh->GetElementGeometry(0) == mfem::Geometry::CUBE, "Expected HEX mesh geometry");

    mfem::BilinearForm mass_form(&m_fes);
    auto *mass_integrator = new mfem::MassIntegrator();

    mfem::IntegrationRule integration_rule = make_hex_tensor_gll_rule_(m_order);
    mass_integrator->SetIntRule(&integration_rule);

    mass_form.AddDomainIntegrator(mass_integrator);
    mass_form.Assemble();
    mass_form.Finalize();

    mfem::SparseMatrix mass_matrix;
    mfem::Array<int> essential_true_dofs;
    mass_form.FormSystemMatrix(essential_true_dofs, mass_matrix);

#ifdef TEST_COMPILE
    std::cout << "\n===== Overlap Matrix =====" << std::endl;
    std::cout << "Size: " << mass_matrix.Height() << " x " << mass_matrix.Width() << std::endl;
    std::cout << "Non-zeros: " << mass_matrix.NumNonZeroElems() << std::endl;

    double off_diagonal_norm = 0.0;
    int off_diagonal_count = 0;
    int negative_diagonal_count = 0;
    const int *row_offsets = mass_matrix.GetI();
    const int *column_indices = mass_matrix.GetJ();
    const double *values = mass_matrix.GetData();
    for (int row = 0; row < mass_matrix.Height(); ++row)
    {
        for (int index = row_offsets[row]; index < row_offsets[row + 1]; ++index)
        {
            const int col = column_indices[index];
            const double value = values[index];
            if (row != col)
            {
                off_diagonal_norm += value * value;
                ++off_diagonal_count;
            }
            else if (value < 0.0)
            {
                ++negative_diagonal_count;
            }
        }
    }

    std::cout << "Off-diagonal elements: " << off_diagonal_count << std::endl;
    std::cout << "Off-diagonal norm: " << std::sqrt(off_diagonal_norm) << std::endl;
    std::cout << "Negative diagonal elements: " << negative_diagonal_count << std::endl;
    std::cout << "Is diagonal: " << (off_diagonal_norm < 1e-10 ? "YES" : "NO") << std::endl;
    std::cout << "\n===== Integration Rule Info =====" << std::endl;
    std::cout << "Order: " << m_order << std::endl;
    std::cout << "Expected GLL points per dimension: " << (m_order + 1) << std::endl;
    std::cout << "Actual integration points: " << integration_rule.GetNPoints() << std::endl;
    std::cout << "\n===== Mesh Info =====" << std::endl;
    std::cout << "Elements: " << m_mesh->GetNE() << std::endl;
    std::cout << "Vertices: " << m_mesh->GetNV() << std::endl;
    std::cout << "DOFs (before periodicity): " << m_fes.GetVSize() << std::endl;
    std::cout << "DOFs (after periodicity): " << m_fes.GetTrueVSize() << std::endl;
#endif

    m_mass_diag_true.SetSize(mass_matrix.Height());
    mass_matrix.GetDiag(m_mass_diag_true);

    m_minv_half_diag_true.SetSize(m_mass_diag_true.Size());
    for (int dof_index = 0; dof_index < m_mass_diag_true.Size(); ++dof_index)
    {
        MFEM_VERIFY(m_mass_diag_true(dof_index) > 0.0, "Mass diagonal must be positive");
        m_minv_half_diag_true(dof_index) = 1.0 / std::sqrt(m_mass_diag_true(dof_index));
    }
}
