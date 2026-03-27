#pragma once
// Keep ndarray on the narrow local-array path:
// - no slicing/subviews
// - use MFEM types instead for MPI-facing or distributed data
#include "ndarray/ndarray.hpp"
#include <complex>

typedef NDArray<double, 1, int> DoubleArray1D;
typedef NDArray<double, 2, int> DoubleArray2D;
typedef NDArray<std::complex<double>, 1, int> ComplexArray1D;
