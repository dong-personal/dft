#pragma once
// do not use slice of ndarray
#include "ndarray/ndarray.hpp"
#include <complex>

typedef NDArray<double, 1, int> DoubleArray1D;
typedef NDArray<double, 2, int> DoubleArray2D;
typedef NDArray<std::complex<double>, 1, int> ComplexArray1D;