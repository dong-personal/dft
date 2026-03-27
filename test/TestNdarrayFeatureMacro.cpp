#include <iostream>

int main()
{
#ifdef __cpp_multidimensional_subscript
    std::cout << "__cpp_multidimensional_subscript = " << __cpp_multidimensional_subscript
              << std::endl;
#if __cpp_multidimensional_subscript >= 202110L
    std::cout << "Supports C++23 multidimensional subscript" << std::endl;
#else
    std::cout << "Does not support C++23 multidimensional subscript" << std::endl;
#endif
#else
    std::cout << "__cpp_multidimensional_subscript is undefined" << std::endl;
#endif

    return 0;
}
