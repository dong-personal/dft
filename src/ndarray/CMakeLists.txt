cmake_minimum_required(VERSION 3.10)
project(test_ndarray)
include(CheckCXXSourceCompiles)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


check_cxx_source_compiles("
int main() {
#if !defined(__cpp_multidimensional_subscript) || __cpp_multidimensional_subscript < 202110L
#error Not supported
#endif
return 0;
}
" HAS_MD_SUBSCRIPT)

if(HAS_MD_SUBSCRIPT)
    message(STATUS "编译器支持 C++23 多维下标运算符")
else()
    message(STATUS "编译器不支持 C++23 多维下标运算符")
endif()

# 查找gtest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# 添加你的源码
add_executable(test_ndarray test_ndarray.cpp)

# 链接gtest和pthread
target_link_libraries(test_ndarray PRIVATE GTest::GTest GTest::Main pthread)