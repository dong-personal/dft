#include "ndarray/ndarray.hpp"

#ifndef TEST_COMPILE
template class NDArray<double, 1, int, 3>;
template class NDArray<double, 2, int, 9>;
#endif
