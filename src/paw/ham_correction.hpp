#pragma once

#ifndef DFT_PAW_HAM_CORRECTION_HPP
#define DFT_PAW_HAM_CORRECTION_HPP
#include "paw/paw_setup.hpp"

namespace dft
{

class HamCorrection
{
  public:
    DenseMatrix &kinetic_energy_differences() { return m_kinetic_energy_differences; }
    const DenseMatrix &kinetic_energy_differences() const { return m_kinetic_energy_differences; }

    DenseMatrix &static_coulomb_correction() { return m_static_coulomb_correction; }
    const DenseMatrix &static_coulomb_correction() const { return m_static_coulomb_correction; }

    DenseMatrix &fixed_nonlocal_correction() { return m_fixed_nonlocal_correction; }
    const DenseMatrix &fixed_nonlocal_correction() const { return m_fixed_nonlocal_correction; }

    void validate(const PAWSetup &setup) const;

  private:
    DenseMatrix m_kinetic_energy_differences;
    DenseMatrix m_static_coulomb_correction;
    DenseMatrix m_fixed_nonlocal_correction;
};

HamCorrection build_ham_correction(const PAWSetup &setup);

} // namespace dft
#endif // DFT_PAW_HAM_CORRECTION_HPP
