add_library(dft_core STATIC
    DFTMesh.cpp
    FEspace.cpp
    PAW.cpp
    kineticOperator.cpp
    pawBasisEvaluator.cpp
    pawNonlocalFixedOperator.cpp
    periodicKDTree.cpp
    radialInterpolator.cpp
    sphericalHarmonics.cpp
    readPoscar.cpp
    ${CMAKE_SOURCE_DIR}/third_party/tinyxml2/tinyxml2.cpp
)

target_include_directories(dft_core
    PUBLIC
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/third_party
        ${CMAKE_SOURCE_DIR}/third_party/tinyxml2
        ${MFEM_INCLUDE_DIR}
)

target_link_libraries(dft_core
    PUBLIC
        ${MFEM_LIBRARY}
)

target_precompile_headers(dft_core PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/pch.hpp"
)
