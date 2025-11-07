cmake_minimum_required(VERSION 3.16)
project(cardillo LANGUAGES C CXX)

add_compile_options(-w)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optional: gprof profiling support
option(ENABLE_GPROF "Enable gprof instrumentation (-pg)" OFF)
if(ENABLE_GPROF)
    message(STATUS "Building with gprof profiling (-pg)")
    add_compile_options(-pg)
    add_link_options(-pg)
endif()

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

# Set release optimization flags
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    message(STATUS "Building with Release optimizations")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        add_compile_options(-O3 -march=native -funroll-loops)
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        add_compile_options(/O2)
    endif()
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

# Prefer system Eigen; fall back to FetchContent when network is available
find_package(Eigen3 3.4...5.0.0 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        Eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 5.0.0
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(Eigen)
else()
    message(STATUS "Found system Eigen3: ${EIGEN3_VERSION}")
    add_compile_definitions(EIGEN3_FOUND)
endif()

# EnTT (header-only) via FetchContent
FetchContent_Declare(
    EnTT
    GIT_REPOSITORY https://github.com/skypjack/entt.git
    GIT_TAG v3.13.2
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(EnTT)

# MPI
find_package(MPI REQUIRED)

# PETSc via pkg-config (optional)
find_package(PkgConfig REQUIRED)
pkg_search_module(PETSC QUIET IMPORTED_TARGET PETSc>=3.15 petsc>=3.15)
if(NOT PETSC_FOUND)
    message(WARNING "PETSc not found. Building without PETSc integration. If you need PETSc, install it or set PKG_CONFIG_PATH so pkg-config can find PETSc.")
endif()

# # Try to find locally-built COAL via pkg-config (coal.pc is provided in coal's build/)
# pkg_search_module(COAL QUIET IMPORTED_TARGET coal)
# if(COAL_FOUND)
#     message(STATUS "Found COAL: ${COAL_VERSION}")
#     add_definitions(-DHAVE_COAL)
#     include_directories(${COAL_INCLUDEDIR})
#     link_directories(${COAL_LIBDIR})
# else()
#     message(WARNING "COAL not found via pkg-config. If you built COAL locally, add its build directory to CMAKE_PREFIX_PATH or PKG_CONFIG_PATH (e.g. -DCMAKE_PREFIX_PATH=/home/$USER/BA/lib/coal-3.0.2/build)")
# endif()

# Python3
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
include_directories( ${PYTHON3_INCLUDE_DIRS})

include(ExternalProject)

# Coal via pip (via ExternalProject_Add)
set(COAL_INSTALL_DIR ${CMAKE_BINARY_DIR}/external/coal)

ExternalProject_Add(coal_project
    PREFIX ${CMAKE_BINARY_DIR}/_deps/coal
    DOWNLOAD_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ${Python3_EXECUTABLE} -m pip install coal -t ${COAL_INSTALL_DIR}
    INSTALL_COMMAND ""
)

# If PETSc was discovered via pkg-config, expose a macro and include dirs
if(PETSC_FOUND)
    message(STATUS "Found PETSc via pkg-config; enabling HAVE_PETSC")
    add_compile_definitions(HAVE_PETSC)
    include_directories(${PETSC_INCLUDE_DIRS})
    link_directories(${PETSC_LIBRARY_DIRS})
endif()

# Add subdirs
add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(benchmarks/matrix_bench)
