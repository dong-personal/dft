add_library(dft_core STATIC
    DFTMesh.cpp
    FEspace.cpp
    kineticOperator.cpp
    readPoscar.cpp
)

target_include_directories(dft_core
    PUBLIC
        ${CMAKE_SOURCE_DIR}/src
        ${MFEM_INCLUDE_DIR}
)

target_link_libraries(dft_core
    PUBLIC
        ${MFEM_LIBRARY}
)

target_precompile_headers(dft_core PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/pch.hpp"
)
