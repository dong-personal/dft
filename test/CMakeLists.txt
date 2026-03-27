set(TEST_SOURCES!
    TestReadPoscar.cpp
    TestDFTMesh.cpp
    TestOverlap.cpp
)

# # -------- test_mesh --------
# add_executable(test_mesh test_mesh.cpp)

# target_include_directories(test_mesh
#     PRIVATE
#         ${CMAKE_SOURCE_DIR}/src
#         ${MFEM_INCLUDE_DIR}
# )

# target_link_libraries(test_mesh
#     PRIVATE
#         dft_core
#         ${MFEM_LIBRARY}
# )

# # -------- test_dft_lda --------
# add_executable(test_lda test_lda.cpp)

# target_include_directories(test_lda
#     PRIVATE
#         ${MFEM_INCLUDE_DIR}
# )

# target_link_libraries(test_lda
#     PRIVATE
#         ${MFEM_LIBRARY}
# )

# # -------- TestReadPoscar --------

# add_executable(TestReadPoscar TestReadPoscar.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../src/readPoscar.cpp)
# target_include_directories(TestReadPoscar
#     PRIVATE
#         ${CMAKE_SOURCE_DIR}/src
# )


# # -------- TestDFTMesh --------
# add_executable(TestDFTMesh TestDFTMesh.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../src/readPoscar.cpp)
# target_include_directories(TestDFTMesh
#     PRIVATE
#         ${CMAKE_SOURCE_DIR}/src
#         ${MFEM_INCLUDE_DIR}
# )
# target_link_libraries(TestDFTMesh
#     PRIVATE
#         ${MFEM_LIBRARY}
# )

# -------- TestOverlap --------
add_executable(TestOverlap TestOverlap.cpp)
target_include_directories(TestOverlap
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${MFEM_INCLUDE_DIR}
)
target_link_libraries(TestOverlap
    PRIVATE
        dft_core
        ${MFEM_LIBRARY}
)
target_compile_definitions(TestOverlap PRIVATE TEST_COMPILE)
target_compile_definitions(TestOverlap PRIVATE DFT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

# -------- TestKinetic --------
add_executable(TestKinetic TestKinetic.cpp)
target_include_directories(TestKinetic
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${MFEM_INCLUDE_DIR}
)
target_link_libraries(TestKinetic
    PRIVATE
        dft_core
        ${MFEM_LIBRARY}
)
target_compile_definitions(TestKinetic PRIVATE TEST_COMPILE)
target_compile_definitions(TestKinetic PRIVATE DFT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

# -------- TestPAW --------
add_executable(TestPAW TestPAW.cpp)
target_include_directories(TestPAW
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${MFEM_INCLUDE_DIR}
)
target_link_libraries(TestPAW
    PRIVATE
        dft_core
        ${MFEM_LIBRARY}
)
target_compile_definitions(TestPAW PRIVATE DFT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

# -------- TestRadialInterpolator --------
add_executable(TestRadialInterpolator TestRadialInterpolator.cpp)
target_include_directories(TestRadialInterpolator
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${MFEM_INCLUDE_DIR}
)
target_link_libraries(TestRadialInterpolator
    PRIVATE
        dft_core
        ${MFEM_LIBRARY}
)
target_compile_definitions(TestRadialInterpolator PRIVATE DFT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

# -------- TestSphericalHarmonics --------
add_executable(TestSphericalHarmonics TestSphericalHarmonics.cpp)
target_include_directories(TestSphericalHarmonics
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${MFEM_INCLUDE_DIR}
)
target_link_libraries(TestSphericalHarmonics
    PRIVATE
        dft_core
        ${MFEM_LIBRARY}
)

# -------- TestPAWBasisEvaluator --------
add_executable(TestPAWBasisEvaluator TestPAWBasisEvaluator.cpp)
target_include_directories(TestPAWBasisEvaluator
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${MFEM_INCLUDE_DIR}
)
target_link_libraries(TestPAWBasisEvaluator
    PRIVATE
        dft_core
        ${MFEM_LIBRARY}
)

# -------- TestPeriodicKDTree --------
add_executable(TestPeriodicKDTree TestPeriodicKDTree.cpp)
target_include_directories(TestPeriodicKDTree
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${MFEM_INCLUDE_DIR}
)
target_link_libraries(TestPeriodicKDTree
    PRIVATE
        dft_core
        ${MFEM_LIBRARY}
)
target_compile_definitions(TestPeriodicKDTree PRIVATE DFT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

# -------- TestPAWNonlocalFixed --------
add_executable(TestPAWNonlocalFixed TestPAWNonlocalFixed.cpp)
target_include_directories(TestPAWNonlocalFixed
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${MFEM_INCLUDE_DIR}
)
target_link_libraries(TestPAWNonlocalFixed
    PRIVATE
        dft_core
        ${MFEM_LIBRARY}
)
target_compile_definitions(TestPAWNonlocalFixed PRIVATE DFT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

# -------- TestNdarray --------
add_executable(TestNdarray TestNdarray.cpp)
target_include_directories(TestNdarray
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
)
target_compile_definitions(TestNdarray PRIVATE TEST_COMPILE)

# -------- TestNdarrayFeatureMacro --------
add_executable(TestNdarrayFeatureMacro TestNdarrayFeatureMacro.cpp)
