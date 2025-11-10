add_library(${PROJECT_NAME} STATIC
    cardillo.cpp
    physics/physics_system.cpp
    physics/assets.cpp
    physics/constraints.cpp
    io/vtk_writer_binary.cpp
    io/vtk_sphere_util.cpp
    io/heightmap_loader.cpp
    io/softbody_loader.cpp
    solver/moreau.cpp
    solver/CPG.cpp
    solver/projected_jacobi.cpp
    misc/spline.cpp
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

if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME}
        PUBLIC
            OpenMP::OpenMP_CXX
    )
endif()

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
