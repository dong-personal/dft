#include <iostream>
int main()
{
#ifdef __cpp_multidimensional_subscript
    std::cout << "__cpp_multidimensional_subscript = " << __cpp_multidimensional_subscript << std::endl;
#if __cpp_multidimensional_subscript >= 202110L
    std::cout << "支持 C++23 多维下标运算符" << std::endl;
#else
    std::cout << "不支持 C++23 多维下标运算符" << std::endl;
#endif
#else
    std::cout << "未定义 __cpp_multidimensional_subscript, 不支持多维下标运算符" << std::endl;
#endif
    return 0;
}