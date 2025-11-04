add_library(${PROJECT_NAME} STATIC
    cardillo.cpp
    physics/physics_system.cpp
    io/vtk_writer_binary.cpp
    io/vtk_sphere_util.cpp
    io/heightmap_loader.cpp
    io/softbody_loader.cpp
    solver/moreau.cpp
    solver/CPG.cpp
    solver/projected_jacobi.cpp
    physics/dynamics_assembler.cpp
    collision/collision_coal.cpp
    partitioning/naive_partitioner.cpp
    comm/communication.cpp
    config/config.cpp
)

target_include_directories(${PROJECT_NAME} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${entt_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/lib/third_party/tinyexr
)

target_link_libraries(${PROJECT_NAME}
    PUBLIC
        Eigen3::Eigen
        PkgConfig::PETSC
        MPI::MPI_CXX
)

# # COAL collision library is required
# # set(COAL_HINT_ROOT "${CMAKE_SOURCE_DIR}/../lib/coal-3.0.2/install")
# set(COAL_HINT_ROOT "/home/jonas/gitprojects/coal_new/")
# if (EXISTS "${COAL_HINT_ROOT}")
#     set(CMAKE_PREFIX_PATH "${COAL_HINT_ROOT};${CMAKE_PREFIX_PATH}")
# endif()
# find_path(COAL_INCLUDE_DIR NAMES coal/collision.h
#     HINTS ${COAL_HINT_ROOT}/include
#     PATH_SUFFIXES include
# )
# # if (NOT COAL_INCLUDE_DIR AND EXISTS "${CMAKE_SOURCE_DIR}/../lib/coal-3.0.2/include")
# #     set(COAL_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/../lib/coal-3.0.2/include")
# # endif()
# find_library(COAL_LIBRARY NAMES coal
#     HINTS ${COAL_HINT_ROOT}/build/lib ${COAL_HINT_ROOT}/build/lib64 #"${CMAKE_SOURCE_DIR}/../lib/coal-3.0.2/build/lib"
# )
# if (NOT COAL_INCLUDE_DIR OR NOT COAL_LIBRARY)
#     message(FATAL_ERROR "COAL is required but not found. Set COAL_INCLUDE_DIR and COAL_LIBRARY or install COAL under lib/coal-3.0.2.")
# endif()
# target_include_directories(${PROJECT_NAME} PUBLIC ${COAL_INCLUDE_DIR})
# # COAL build tree may also generate headers in build/include; add it if present
# # if (EXISTS "${CMAKE_SOURCE_DIR}/../lib/coal-3.0.2/build/include")
# #     target_include_directories(${PROJECT_NAME} PUBLIC "${CMAKE_SOURCE_DIR}/../lib/coal-3.0.2/build/include")
# # endif()
# if (EXISTS ${COAL_HINT_ROOT}/build/include)
#     target_include_directories(${PROJECT_NAME} PUBLIC ${COAL_HINT_ROOT}/build/include)
# endif()
# target_link_libraries(${PROJECT_NAME} PUBLIC ${COAL_LIBRARY})
# target_compile_definitions(${PROJECT_NAME} PUBLIC CARDILLO_HAVE_COAL=1)

add_dependencies(${PROJECT_NAME} coal_project)

# Include and link directories after install
target_include_directories(${PROJECT_NAME} PUBLIC
    "${COAL_INSTALL_DIR}/cmeel.prefix/include/"
)

target_link_directories(${PROJECT_NAME} PUBLIC
    "${COAL_INSTALL_DIR}/cmeel.prefix/lib/"
)

target_link_libraries(${PROJECT_NAME} INTERFACE coal)

target_compile_definitions(${PROJECT_NAME} PUBLIC CARDILLO_HAVE_COAL=1)

# Optional: enable EXR loading via tinyexr (header-only). Define CARDILLO_USE_TINYEXR to turn it on if header is provided.
# Users can set -DCARDILLO_USE_TINYEXR=ON when configuring to enable.
option(CARDILLO_USE_TINYEXR "Enable EXR heightmap loading via tinyexr" OFF)
if (CARDILLO_USE_TINYEXR)
    target_compile_definitions(${PROJECT_NAME} PUBLIC CARDILLO_USE_TINYEXR=1)
    find_package(ZLIB REQUIRED)
    target_link_libraries(${PROJECT_NAME} PUBLIC ZLIB::ZLIB)
endif()

# Simple tests (no framework): build small executables
add_executable(matrix_tests tests/matrix_tests.cpp)
target_link_libraries(matrix_tests PRIVATE ${PROJECT_NAME})
