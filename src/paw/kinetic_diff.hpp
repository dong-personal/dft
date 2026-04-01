#pragma once
#ifndef DFT_PAW_KINETIC_DIFF_HPP
#include "paw/paw_setup.hpp"

#include <cstddef>
#include <vector>

namespace dft
{

DenseMatrix build_full_kinetic_diff_matrix(const std::vector<double> &flat_values, std::size_t state_count);

} // namespace dft
#endif