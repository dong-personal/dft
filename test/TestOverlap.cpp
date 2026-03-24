#include "FEspace.h"
#include <cmath>
#include <iostream>
#include <memory>

int main(void)
{
    const std::string poscar_path = std::string(DFT_SOURCE_DIR) + "/data/C.POSCAR";

    auto structure = std::make_shared<dft::Structure>(dft::readPoscar(poscar_path));
    auto dft_mesh = std::make_shared<dft::DFTMesh>();
    dft_mesh->set_structure(structure);
    dft_mesh->init_periodic_cell_from_lattice(structure->lattice(), 10, 10, 10);
    // dft_mesh->refine_near_atoms(1.0, 2, 0.5);

    DFTGLLHexSpace fespace(dft_mesh, 3);

    const mfem::Vector &diag = fespace.MassDiagTrue();
    const mfem::Vector &minvhalf = fespace.MinvHalfDiagTrue();

    if (diag.Size() != fespace.getDOF())
    {
        std::cerr << "Mass diagonal size mismatch: " << diag.Size() << " vs true dofs " << fespace.getDOF()
                  << std::endl;
        return 1;
    }

    if (minvhalf.Size() != diag.Size())
    {
        std::cerr << "M^{-1/2} size mismatch: " << minvhalf.Size() << " vs " << diag.Size() << std::endl;
        return 1;
    }

    double min_diag = diag(0);
    double max_diag = diag(0);
    double identity_error = 0.0;
    for (int i = 0; i < diag.Size(); ++i)
    {
        if (diag(i) <= 0.0)
        {
            std::cerr << "Non-positive mass diagonal entry at " << i << ": " << diag(i) << std::endl;
            return 1;
        }

        min_diag = std::min(min_diag, diag(i));
        max_diag = std::max(max_diag, diag(i));

        const double recovered = diag(i) * minvhalf(i) * minvhalf(i);
        const double err = std::abs(recovered - 1.0);
        identity_error = std::max(identity_error, err);
    }

    mfem::Vector ones(diag.Size());
    ones = 1.0;

    mfem::Vector mass_times_ones;
    mfem::Vector scaled_ones;
    fespace.ApplyMassDiagTrue(ones, mass_times_ones);
    fespace.ApplyMinvHalfTrue(ones, scaled_ones);

    std::cout << "Verified overlap diagonal on true dofs" << std::endl;
    std::cout << "True dofs: " << diag.Size() << std::endl;
    std::cout << "Min mass entry: " << min_diag << std::endl;
    std::cout << "Max mass entry: " << max_diag << std::endl;
    std::cout << "Max error in M * (M^{-1/2})^2 = I: " << identity_error << std::endl;
    std::cout << "||M * 1||_2 = " << mass_times_ones.Norml2() << std::endl;
    std::cout << "||M^{-1/2} * 1||_2 = " << scaled_ones.Norml2() << std::endl;

    if (identity_error > 1e-12)
    {
        std::cerr << "Inverse square root mass diagonal failed consistency check" << std::endl;
        return 1;
    }

    return 0;
}
