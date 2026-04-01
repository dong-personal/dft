#pragma once
#ifndef DFT_PAW_COULOMB_CORRECTION_HPP
#include "paw/paw_setup.hpp"

namespace dft
{

DenseMatrix build_two_index_coulomb_correction(const PAWSetup &setup);

} // namespace dft

#endif
