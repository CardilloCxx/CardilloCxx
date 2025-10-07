add_library(${PROJECT_NAME} STATIC mylib.cpp)

target_include_directories(${PROJECT_NAME} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(mylib
    PUBLIC
        Eigen3::Eigen
        PkgConfig::PETSC
        MPI::MPI_CXX
)
