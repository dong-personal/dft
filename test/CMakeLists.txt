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
add_executable(TestOverlap TestOverlap.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../src/readPoscar.cpp)
target_include_directories(TestOverlap
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${MFEM_INCLUDE_DIR}
)
target_link_libraries(TestOverlap
    PRIVATE
        ${MFEM_LIBRARY}
)
target_compile_definitions(TestOverlap PRIVATE TEST_COMPILE)