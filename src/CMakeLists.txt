add_library(${PROJECT_NAME} STATIC
    cardillo.cpp
    physics/physics_system.cpp
    io/vtk_writer.cpp
    solver/moreau.cpp
    solver/CPG.cpp
    solver/projected_jacobi.cpp
    physics/dynamics_assembler.cpp
    collision/broad_phase.cpp
    collision/narrow_phase.cpp
    collision/box_box_collision.cpp
    collision/collision_manager.cpp
    partitioning/naive_partitioner.cpp
    comm/communication.cpp
)

target_include_directories(${PROJECT_NAME} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${entt_SOURCE_DIR}/src
)

target_link_libraries(${PROJECT_NAME}
    PUBLIC
        Eigen3::Eigen
        PkgConfig::PETSC
        MPI::MPI_CXX
)

# Simple tests (no framework): build small executables
add_executable(matrix_tests tests/matrix_tests.cpp)
target_link_libraries(matrix_tests PRIVATE ${PROJECT_NAME})
