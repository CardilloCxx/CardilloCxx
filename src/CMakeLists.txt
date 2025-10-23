add_library(${PROJECT_NAME} STATIC
    cardillo.cpp
    physics/physics_system.cpp
    io/vtk_writer.cpp
    solver/moreau.cpp
    solver/CPG.cpp
    solver/projected_jacobi.cpp
    physics/dynamics_assembler.cpp
    # deprecated collision pipeline moved under collision/depracted (not built)
    collision/collision_coal.cpp
    partitioning/naive_partitioner.cpp
    comm/communication.cpp
    config/config.cpp
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

# COAL collision library is required
set(COAL_HINT_ROOT "${CMAKE_SOURCE_DIR}/../lib/coal-3.0.2/install")
if (EXISTS "${COAL_HINT_ROOT}")
    set(CMAKE_PREFIX_PATH "${COAL_HINT_ROOT};${CMAKE_PREFIX_PATH}")
endif()
find_path(COAL_INCLUDE_DIR NAMES coal/collision.h
    HINTS ${COAL_HINT_ROOT}/include
    PATH_SUFFIXES include
)
if (NOT COAL_INCLUDE_DIR AND EXISTS "${CMAKE_SOURCE_DIR}/../lib/coal-3.0.2/include")
    set(COAL_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/../lib/coal-3.0.2/include")
endif()
find_library(COAL_LIBRARY NAMES coal
    HINTS ${COAL_HINT_ROOT}/lib ${COAL_HINT_ROOT}/lib64 "${CMAKE_SOURCE_DIR}/../lib/coal-3.0.2/build/lib"
)
if (NOT COAL_INCLUDE_DIR OR NOT COAL_LIBRARY)
    message(FATAL_ERROR "COAL is required but not found. Set COAL_INCLUDE_DIR and COAL_LIBRARY or install COAL under lib/coal-3.0.2.")
endif()
target_include_directories(${PROJECT_NAME} PUBLIC ${COAL_INCLUDE_DIR})
# COAL build tree may also generate headers in build/include; add it if present
if (EXISTS "${CMAKE_SOURCE_DIR}/../lib/coal-3.0.2/build/include")
    target_include_directories(${PROJECT_NAME} PUBLIC "${CMAKE_SOURCE_DIR}/../lib/coal-3.0.2/build/include")
endif()
target_link_libraries(${PROJECT_NAME} PUBLIC ${COAL_LIBRARY})
target_compile_definitions(${PROJECT_NAME} PUBLIC CARDILLO_HAVE_COAL=1)

# Simple tests (no framework): build small executables
add_executable(matrix_tests tests/matrix_tests.cpp)
target_link_libraries(matrix_tests PRIVATE ${PROJECT_NAME})
