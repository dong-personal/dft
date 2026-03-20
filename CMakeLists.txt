cmake_minimum_required(VERSION 3.12)
project(DFT_FEM_Project CXX)
include(CheckCXXSourceCompiles)

set(CMAKE_CXX_STANDARD 23) 
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


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

# ------------ 手动指定 MFEM 根目录、include 和 lib 目录 ------------

set(MFEM_ROOT "${CMAKE_SOURCE_DIR}/mfem")

# 头文件目录：mfem/include
set(MFEM_INCLUDE_DIR "${MFEM_ROOT}/include")

# 库文件：根据你实际的文件名来改！
# 如果是 Linux/macOS 通常是 libmfem.a 或 libmfem.so
set(MFEM_LIBRARY "${MFEM_ROOT}/lib/libmfem.a")

# 如果是 Windows，可能是这样：
# set(MFEM_LIBRARY "${MFEM_ROOT}/lib/mfem.lib")

# 检查一下路径是否存在，方便尽早报错
if (NOT EXISTS "${MFEM_LIBRARY}")
    message(FATAL_ERROR "MFEM library not found at: ${MFEM_LIBRARY}")
endif()

if (NOT EXISTS "${MFEM_INCLUDE_DIR}/mfem.hpp")
    message(FATAL_ERROR "mfem.hpp not found in: ${MFEM_INCLUDE_DIR}")
endif()

message(STATUS "Using MFEM include: ${MFEM_INCLUDE_DIR}")
message(STATUS "Using MFEM library: ${MFEM_LIBRARY}")

# 让后面的 target 可以用这两个变量
set(MFEM_INCLUDE_DIR ${MFEM_INCLUDE_DIR} CACHE PATH "MFEM include dir")
set(MFEM_LIBRARY     ${MFEM_LIBRARY}     CACHE FILEPATH "MFEM library")

add_subdirectory(src)
add_subdirectory(test)
