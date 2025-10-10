add_library(${PROJECT_NAME} STATIC cardillo.cpp physics/physics_system.cpp io/vtk_writer.cpp solver/moreau.cpp)

target_include_directories(${PROJECT_NAME} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(${PROJECT_NAME}
    PUBLIC
        Eigen3::Eigen
        PkgConfig::PETSC
        MPI::MPI_CXX
)
