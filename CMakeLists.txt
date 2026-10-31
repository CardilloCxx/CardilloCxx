# 3.24: FetchContent with FIND_PACKAGE_ARGS (cmake/Dependencies.cmake), header file sets.
cmake_minimum_required(VERSION 3.24)
project(cardillo VERSION 0.1.0 LANGUAGES C CXX)

if(NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 20)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

if(PROJECT_IS_TOP_LEVEL)
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
    # Without a build type CMake compiles without optimization; default to RelWithDebInfo.
    if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
        set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Build type" FORCE)
        set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
    endif()
endif()

# --- Options ---
# Examples, tests, install rules and -march=native default to ON only when CardilloCxx is the
# top-level project, so that a parent project using it through add_subdirectory/FetchContent gets
# just the library.
option(CARDILLO_BUILD_EXAMPLES "Build the example scene executables" ${PROJECT_IS_TOP_LEVEL})
option(CARDILLO_BUILD_TESTS "Build the tests (ctest)" ${PROJECT_IS_TOP_LEVEL})
option(CARDILLO_INSTALL "Generate install rules and the CMake package (find_package(cardillo))" ${PROJECT_IS_TOP_LEVEL})
# -march=native ties the binaries to the build machine's CPU; turn it off for binaries run elsewhere.
option(CARDILLO_ENABLE_NATIVE_ARCH "Compile with -march=native" ${PROJECT_IS_TOP_LEVEL})
option(CARDILLO_ENABLE_TSAN "Build with ThreadSanitizer (use with CMAKE_BUILD_TYPE=Debug)" OFF)
option(CARDILLO_ENABLE_ASAN "Build with AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
# Use an installed Coal instead of fetching the fork; it must be built from the fork (custom shapes).
option(CARDILLO_USE_SYSTEM_COAL "Use an installed Coal built from the 41pha1/coal fork" OFF)

# Optional conic solver backends (solver.type = clarabel | qoco | conicxx | scs | mosek). The production
# solvers (projected Jacobi, condensed) and the reference solvers are always built.
option(CARDILLO_WITH_CLARABEL "Clarabel interior-point backend (needs Rust/cargo)" OFF)
option(CARDILLO_WITH_QOCO "QOCO interior-point backend" OFF)
option(CARDILLO_WITH_CONICXX "ConicXX interior-point backend" OFF)
option(CARDILLO_WITH_SCS "SCS ADMM backend" OFF)
# MOSEK is not fetched: it needs an installation and a license (cmake/FindMOSEK.cmake, MOSEK_ROOT).
option(CARDILLO_WITH_MOSEK "MOSEK interior-point backend (needs a MOSEK installation and license)" OFF)
option(CARDILLO_SCS_OPENMP "Build SCS with its OpenMP-parallel cone projections" ON)

# --- Compile flags of this project and of the dependencies it fetches ---
if(CARDILLO_ENABLE_TSAN AND CARDILLO_ENABLE_ASAN)
    message(FATAL_ERROR "CARDILLO_ENABLE_TSAN and CARDILLO_ENABLE_ASAN cannot be combined.")
endif()
if(CARDILLO_ENABLE_TSAN)
    add_compile_options(-fsanitize=thread -g -O1)
    add_link_options(-fsanitize=thread)
endif()
if(CARDILLO_ENABLE_ASAN)
    add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer -g)
    add_link_options(-fsanitize=address,undefined)
endif()
if(CARDILLO_ENABLE_NATIVE_ARCH)
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag("-march=native" CARDILLO_HAS_MARCH_NATIVE)
    if(CARDILLO_HAS_MARCH_NATIVE)
        # Global, so that the fetched dependencies use the same ISA; it is also a usage requirement of
        # cardillo (src/CMakeLists.txt) because it changes the alignment of Eigen types in the API.
        add_compile_options(-march=native)
    else()
        message(WARNING "The compiler does not support -march=native; building for the default ISA.")
    endif()
endif()

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include(cmake/Dependencies.cmake)
include(cmake/ConicBackends.cmake)

add_subdirectory(src)
if(CARDILLO_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()
if(CARDILLO_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
if(CARDILLO_INSTALL)
    include(cmake/Install.cmake)
endif()

message(STATUS "CardilloCxx ${PROJECT_VERSION}: Clarabel=${CARDILLO_WITH_CLARABEL} QOCO=${CARDILLO_WITH_QOCO} ConicXX=${CARDILLO_WITH_CONICXX} SCS=${CARDILLO_WITH_SCS} MOSEK=${CARDILLO_WITH_MOSEK}")
