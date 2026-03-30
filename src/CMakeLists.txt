add_library(dft_local_core STATIC
    paw/paw.cpp
    ndarray/ndarray_instantiations.cpp
    paw/paw_basis_evaluator.cpp
    paw/radial_interpolator.cpp
    paw/spherical_harmonics.cpp
    read_poscar.cpp
    ${CMAKE_SOURCE_DIR}/third_party/tinyxml2/tinyxml2.cpp
)

target_include_directories(dft_local_core
    PUBLIC
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/third_party
        ${CMAKE_SOURCE_DIR}/third_party/tinyxml2
)

target_precompile_headers(dft_local_core PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/local_pch.hpp"
)

add_library(dft_mfem_core STATIC
    dftmesh.cpp
    fespace.cpp
    kinetic_operator.cpp
    paw/paw_nonlocal_fixed_operator.cpp
    periodic_kd_tree.cpp
)

target_include_directories(dft_mfem_core
    PUBLIC
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/third_party
        ${CMAKE_SOURCE_DIR}/third_party/tinyxml2
        ${MFEM_INCLUDE_DIR}
)

target_link_libraries(dft_mfem_core
    PUBLIC
        dft_local_core
        ${MFEM_LIBRARY}
)

target_precompile_headers(dft_mfem_core PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/pch.hpp"
)

add_library(dft_core INTERFACE)

target_link_libraries(dft_core
    INTERFACE
        dft_mfem_core
)
