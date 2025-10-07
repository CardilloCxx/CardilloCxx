cmake_minimum_required(VERSION 3.16)
project(mylib LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set a default build type if none was specified
# see https://blog.kitware.com/cmake-and-the-default-build-type/
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
    STRING "Choose the type of build." FORCE)
# Set the possible values of build type for cmake-gui
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build 
        directory) and run CMake from there. You may need to remove CMakeCache.txt.")
endif()

# nicely support folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# enables usage of standard install directories as ${CMAKE_INSTALL_INCLUDEDIR}
include(GNUInstallDirs)

include(FetchContent)

# Eigen via FetchContent
FetchContent_Declare(
    Eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG master
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(Eigen)

# MPI
find_package(MPI REQUIRED)

# PETSc via pkg-config
find_package(PkgConfig REQUIRED)
pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc>=3.15 petsc>=3.15)

# Add subdirs
add_subdirectory(src)
add_subdirectory(examples)
