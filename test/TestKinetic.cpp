#include "fespace.h"
#include "kinetic_operator.hpp"

#include <cmath>
#include <iostream>
#include <memory>

int main(void)
{
    const std::string poscar_path = std::string(DFT_SOURCE_DIR) + "/data/C.POSCAR";

    auto structure = std::make_shared<dft::Structure>(dft::read_poscar(poscar_path));
    auto dft_mesh = std::make_shared<dft::DFTMesh>();
    dft_mesh->set_structure(structure);
    dft_mesh->init_periodic_cell_from_lattice(structure->lattice(), 10, 10, 10);

    DFTGLLHexSpace fespace(dft_mesh, 3);
    KineticOperator kinetic_operator(fespace);
    kinetic_operator.Assemble();

    const mfem::SparseMatrix &kinetic = kinetic_operator.Matrix();
    const int ndof = fespace.getDOF();
    if (kinetic.Height() != ndof || kinetic.Width() != ndof)
    {
        std::cerr << "Kinetic matrix size mismatch: " << kinetic.Height() << " x " << kinetic.Width() << " for "
                  << ndof << " true dofs" << std::endl;
        return 1;
    }

    mfem::Vector ones(ndof);
    ones = 1.0;

    mfem::Vector kinetic_times_ones;
    kinetic_operator.Mult(ones, kinetic_times_ones);

    const double constant_residual = kinetic_times_ones.Norml2();
    const double constant_energy = kinetic_operator.Energy(ones);

    mfem::Vector probe(ndof);
    for (int i = 0; i < ndof; ++i)
    {
        probe(i) = 1.0 + 0.01 * std::sin(0.1 * static_cast<double>(i));
    }

    const double probe_energy = kinetic_operator.Energy(probe);

    std::cout << "Verified kinetic term on true dofs" << std::endl;
    std::cout << "True dofs: " << ndof << std::endl;
    std::cout << "Kinetic matrix nnz: " << kinetic.NumNonZeroElems() << std::endl;
    std::cout << "||T * 1||_2 = " << constant_residual << std::endl;
    std::cout << "<1|T|1> = " << constant_energy << std::endl;
    std::cout << "<probe|T|probe> = " << probe_energy << std::endl;

    if (constant_energy < -1e-10)
    {
        std::cerr << "Kinetic energy should be non-negative, got " << constant_energy << std::endl;
        return 1;
    }

    if (probe_energy < -1e-10)
    {
        std::cerr << "Probe kinetic energy should be non-negative, got " << probe_energy << std::endl;
        return 1;
    }

    if (constant_residual > 1e-8)
    {
        std::cerr << "Constant state should be close to the null space of the periodic kinetic term, got ||T*1|| = "
                  << constant_residual << std::endl;
        return 1;
    }

    return 0;
}
