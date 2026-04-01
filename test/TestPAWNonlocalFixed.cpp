#include "fespace.h"
#include "paw/grid_map.hpp"
#include "paw/ham_correction.hpp"
#include "paw/paw_setup.hpp"
#include "paw/paw_nonlocal_fixed_operator.hpp"

#include <cmath>
#include <iostream>
#include <memory>

namespace
{

constexpr double kMatrixTolerance = 1e-10;

std::size_t magnetic_channel_count(const dft::PAWSetup &setup)
{
    std::size_t count = 0;
    for (const dft::PAWState &state : setup.states())
    {
        count += static_cast<std::size_t>(2 * state.l + 1);
    }

    return count;
}

dft::PeriodicGridPointLocator::Point3 to_locator_point(const dft::Atom::AtomicPosition &point)
{
    dft::PeriodicGridPointLocator::Point3 locator_point(typename dft::PeriodicGridPointLocator::Point3::ShapeType{3});
    locator_point[0] = point[0];
    locator_point[1] = point[1];
    locator_point[2] = point[2];
    return locator_point;
}

void apply_dense_matrix(const dft::DenseMatrix &matrix, const mfem::Vector &x, mfem::Vector &y)
{
    y.SetSize(x.Size());
    y = 0.0;

    for (int row = 0; row < x.Size(); ++row)
    {
        double value = 0.0;
        for (int column = 0; column < x.Size(); ++column)
        {
            value += matrix[row, column] * x(column);
        }

        y(row) = value;
    }
}

double symmetry_residual(const dft::DenseMatrix &matrix)
{
    double residual = 0.0;
    const int extent = matrix.shape()[0];
    for (int row = 0; row < extent; ++row)
    {
        for (int column = 0; column < extent; ++column)
        {
            const double difference = matrix[row, column] - matrix[column, row];
            residual += difference * difference;
        }
    }

    return std::sqrt(residual);
}

} // namespace

int main()
{
    const std::string poscar_path = std::string(DFT_SOURCE_DIR) + "/data/C.POSCAR";
    const std::string paw_path = std::string(DFT_SOURCE_DIR) + "/data/C.GGA_PBE-JTH.xml";

    auto structure = std::make_shared<dft::Structure>(dft::read_poscar(poscar_path));
    auto dft_mesh = std::make_shared<dft::DFTMesh>();
    dft_mesh->set_structure(structure);
    dft_mesh->init_periodic_cell_from_lattice(structure->lattice(), 4, 4, 4);

    DFTGLLHexSpace fespace(dft_mesh, 1);
    const dft::PAWSetup setup = dft::load_paw_setup_xml(paw_path);
    const dft::HamCorrection correction = dft::build_ham_correction(setup);
    const dft::Atom &atom_species = structure->atoms().front();
    const std::size_t position_index = 0;
    const dft::Atom::AtomicPosition atom_center = atom_species.position(position_index);

    PAWNonlocalFixedOperator nonlocal_fixed(fespace, setup, correction, atom_species, position_index);
    dft::PeriodicGridPointLocator locator(fespace);
    const std::vector<dft::PeriodicNeighbor> neighbors =
        locator.radius_search(to_locator_point(atom_center), setup.cutoff_radius());
    const dft::ProjectorTable projector_table =
        dft::build_projector_table(setup, fespace.mesh_source().lattice(), atom_center, locator, neighbors);
    const dft::AtomGridMap grid_map =
        dft::build_atom_grid_map(correction, projector_table, neighbors, fespace.MassDiagTrue());
    dft::PAWSetupRegistry registry;
    registry.add(setup);
    const dft::AtomGLLHamCorrectionMatrix atom_gll_correction(fespace, setup, correction, atom_species, position_index);
    const dft::DenseMatrix &cached_atom_matrix = atom_gll_correction.matrix();
    const dft::DenseMatrix atom_matrix =
        dft::build_atom_gll_ham_correction_matrix(fespace, setup, correction, atom_species, position_index);
    const dft::DenseMatrix total_matrix = dft::build_gll_ham_correction_matrix(fespace, registry);

    mfem::Vector x(fespace.getDOF());
    x = 1.0;
    mfem::Vector y;
    nonlocal_fixed.apply(x, y);
    mfem::Vector y_from_grid_map;
    dft::apply_atom_grid_map(grid_map, x, y_from_grid_map);
    mfem::Vector y_from_atom_matrix;
    apply_dense_matrix(atom_matrix, x, y_from_atom_matrix);
    mfem::Vector y_from_cached_atom_matrix;
    apply_dense_matrix(cached_atom_matrix, x, y_from_cached_atom_matrix);
    mfem::Vector y_from_total_matrix;
    apply_dense_matrix(total_matrix, x, y_from_total_matrix);

    std::cout << "True dofs: " << fespace.getDOF() << '\n';
    std::cout << "Projectors: " << setup.projectors_by_state().size() << '\n';
    std::cout << "Fixed matrix size: " << correction.fixed_nonlocal_correction().shape()[0] << '\n';
    std::cout << "Neighbor hits: " << nonlocal_fixed.neighbors().size() << '\n';
    std::cout << "Mapped grid points: " << grid_map.point_count() << '\n';
    std::cout << "GLL matrix size: " << total_matrix.shape()[0] << '\n';
    std::cout << "||V_nonlocal^fixed * 1||_2 = " << y.Norml2() << '\n';

    if (static_cast<std::size_t>(correction.fixed_nonlocal_correction().shape()[0]) != magnetic_channel_count(setup))
    {
        std::cerr << "Fixed nonlocal matrix size does not match projector count" << std::endl;
        return 1;
    }

    if (nonlocal_fixed.neighbors().empty())
    {
        std::cerr << "No nearby grid points were found for the PAW nonlocal operator" << std::endl;
        return 1;
    }

    if (grid_map.point_count() != neighbors.size())
    {
        std::cerr << "Atom grid map point count does not match neighbor count" << std::endl;
        return 1;
    }

    if (y.Size() != fespace.getDOF())
    {
        std::cerr << "Operator output size mismatch" << std::endl;
        return 1;
    }

    if (y_from_grid_map.Size() != fespace.getDOF())
    {
        std::cerr << "Grid-map output size mismatch" << std::endl;
        return 1;
    }

    if (static_cast<int>(atom_matrix.shape()[0]) != fespace.getDOF() ||
        static_cast<int>(atom_matrix.shape()[1]) != fespace.getDOF() ||
        static_cast<int>(cached_atom_matrix.shape()[0]) != fespace.getDOF() ||
        static_cast<int>(cached_atom_matrix.shape()[1]) != fespace.getDOF() ||
        static_cast<int>(total_matrix.shape()[0]) != fespace.getDOF() ||
        static_cast<int>(total_matrix.shape()[1]) != fespace.getDOF())
    {
        std::cerr << "GLL Hamiltonian correction matrix size mismatch" << std::endl;
        return 1;
    }

    y_from_grid_map -= y;
    if (y_from_grid_map.Norml2() > 1e-12)
    {
        std::cerr << "Grid-map application does not match PAWNonlocalFixedOperator" << std::endl;
        return 1;
    }

    if (symmetry_residual(atom_matrix) > kMatrixTolerance)
    {
        std::cerr << "Atom GLL Hamiltonian correction matrix is not symmetric: "
                  << symmetry_residual(atom_matrix) << std::endl;
        return 1;
    }

    y_from_cached_atom_matrix -= y_from_atom_matrix;
    if (y_from_cached_atom_matrix.Norml2() > 1e-12)
    {
        std::cerr << "Cached atom GLL Hamiltonian correction matrix does not match free-function build" << std::endl;
        return 1;
    }

    if (symmetry_residual(total_matrix) > kMatrixTolerance)
    {
        std::cerr << "Full GLL Hamiltonian correction matrix is not symmetric: "
                  << symmetry_residual(total_matrix) << std::endl;
        return 1;
    }

    if (y_from_atom_matrix.Norml2() == 0.0 || y_from_total_matrix.Norml2() == 0.0)
    {
        std::cerr << "GLL Hamiltonian correction matrix application is identically zero" << std::endl;
        return 1;
    }

    return 0;
}
