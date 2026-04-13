cmake_minimum_required(VERSION 3.16)
project(cardillo LANGUAGES C CXX)

add_compile_options(-w)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY TRUE)

# Optional: gprof profiling support
option(ENABLE_GPROF "Enable gprof instrumentation (-pg)" OFF)
if(ENABLE_GPROF)
    message(STATUS "Building with gprof profiling (-pg)")
    add_compile_options(-pg)
    add_link_options(-pg)
endif()

# set a default build type if none was specified
# see https://blog.kitware.com/cmake-and-the-default-build-type/

# --- OpenMP ---
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found: enabling support")
    add_compile_options(${OpenMP_CXX_FLAGS})
    add_link_options(${OpenMP_CXX_FLAGS})
else()
    message(WARNING "OpenMP not found!")
endif()

set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
    STRING "Choose the type of build." FORCE)
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
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

# Clarabel.rs build toolchain check (used by Clarabel.cpp rust_wrapper).
# Clarabel's CMake uses plain `cargo` in custom commands, so ensure PATH for build tools.
find_program(CARGO_EXECUTABLE cargo HINTS "$ENV{HOME}/.cargo/bin")
if(NOT CARGO_EXECUTABLE)
    message(WARNING
        "Cargo was not found in PATH. Clarabel build requires Rust/Cargo.\n"
        "Install with rustup: curl https://sh.rustup.rs -sSf | sh\n"
        "Then source ~/.cargo/env and reconfigure/build.")
else()
    message(STATUS "Found Cargo: ${CARGO_EXECUTABLE}")

    get_filename_component(CARGO_BIN_DIR "${CARGO_EXECUTABLE}" DIRECTORY)
    set(_cargo_rule_launcher "${CMAKE_COMMAND} -E env \"PATH=${CARGO_BIN_DIR}:$ENV{PATH}\"")
    set_property(GLOBAL PROPERTY RULE_LAUNCH_CUSTOM "${_cargo_rule_launcher}")
endif()

# Build both QOCO backends to allow runtime backend selection.
option(QOCO_USE_CUDA_BACKEND "Build QOCO CUDA backend variant" ON)

set(QOCO_RUNTIME_DIRS "")
set(QOCO_HAVE_CUDA_BACKEND OFF)

if(QOCO_USE_CUDA_BACKEND)
    find_package(CUDAToolkit QUIET)

    find_path(QOCO_CUDSS_INCLUDE cudss.h
        HINTS ${QOCO_CUDSS_INCLUDE_DIR} ${CUDAToolkit_INCLUDE_DIRS}
    )

    find_library(QOCO_CUDSS_LIB
        NAMES cudss
        HINTS ${QOCO_CUDSS_LIBRARY_DIR} ${CUDAToolkit_LIBRARY_DIR}
    )

    if(CUDAToolkit_FOUND AND QOCO_CUDSS_INCLUDE AND QOCO_CUDSS_LIB)
        set(QOCO_HAVE_CUDA_BACKEND ON)

        get_filename_component(QOCO_CUDSS_LIBRARY_DIR "${QOCO_CUDSS_LIB}" DIRECTORY)

        list(APPEND QOCO_RUNTIME_DIRS
            ${QOCO_CUDSS_LIBRARY_DIR}
            ${CUDAToolkit_LIBRARY_DIR}
        )

        list(REMOVE_DUPLICATES QOCO_RUNTIME_DIRS)

        message(STATUS "QOCO CUDA backend variant enabled")
    else()
        message(WARNING "CUDA backend requested but cuDSS not found, CUDA variant will not be built.")
    endif()
endif()

message(STATUS "QOCO CPU backend variant enabled")
if(QOCO_HAVE_CUDA_BACKEND)
    message(STATUS "QOCO CUDA backend variant enabled")
endif()

if(QOCO_ALGEBRA_BACKEND_VALUE STREQUAL "cuda" AND QOCO_RUNTIME_DIRS)
    list(APPEND CMAKE_BUILD_RPATH ${QOCO_RUNTIME_DIRS})
endif()

# Auto-detect QOCO if a local checkout exists at project root
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/qoco")
    message(STATUS "Found local qoco directory, enabling QOCO integration")
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/qoco EXCLUDE_FROM_ALL)
    add_compile_definitions(HAVE_QOCO=1)
endif()

# Coal via pip (via ExternalProject_Add)
set(COAL_INSTALL_DIR ${CMAKE_BINARY_DIR}/external/coal)

ExternalProject_Add(coal_project
    PREFIX ${CMAKE_BINARY_DIR}/_deps/coal
    DOWNLOAD_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ${Python3_EXECUTABLE} -m pip install coal -t ${COAL_INSTALL_DIR}
    INSTALL_COMMAND ""
)

# Install QOCO variants into the active build tree to avoid cross-build contamination.
set(QOCO_INSTALL_DIR_CPU "${CMAKE_BINARY_DIR}/qoco-install-cpu" CACHE PATH "QOCO CPU install prefix")
set(QOCO_INSTALL_DIR_CUDA "${CMAKE_BINARY_DIR}/qoco-install-cuda" CACHE PATH "QOCO CUDA install prefix")

set(QOCO_CMAKE_ARGS_COMMON
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DQOCO_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DBUILD_SHARED_LIBS=ON
    -DENABLE_TESTING=OFF
)

set(QOCO_CMAKE_ARGS_CPU
    ${QOCO_CMAKE_ARGS_COMMON}
    -DCMAKE_INSTALL_PREFIX=${QOCO_INSTALL_DIR_CPU}
    -DQOCO_ALGEBRA_BACKEND=builtin
)

set(QOCO_CMAKE_ARGS_CUDA
    ${QOCO_CMAKE_ARGS_COMMON}
    -DCMAKE_INSTALL_PREFIX=${QOCO_INSTALL_DIR_CUDA}
    -DQOCO_ALGEBRA_BACKEND=cuda
)

if(QOCO_CUDSS_LIBRARY_DIR)
    list(APPEND QOCO_CMAKE_ARGS_CUDA -DCMAKE_LIBRARY_PATH=${QOCO_CUDSS_LIBRARY_DIR})
endif()
if(QOCO_CUDSS_INCLUDE_DIR)
    list(APPEND QOCO_CMAKE_ARGS_CUDA -DCMAKE_INCLUDE_PATH=${QOCO_CUDSS_INCLUDE_DIR})
endif()
if(QOCO_RUNTIME_DIRS)
    string(JOIN ":" QOCO_EXTERNAL_RPATH ${QOCO_RUNTIME_DIRS})
    list(APPEND QOCO_CMAKE_ARGS_CUDA
        -DCMAKE_BUILD_RPATH=${QOCO_EXTERNAL_RPATH}
        -DCMAKE_INSTALL_RPATH=${QOCO_EXTERNAL_RPATH}
    )
endif()

ExternalProject_Add(qoco_project_cpu
    GIT_REPOSITORY https://github.com/qoco-org/qoco.git
    GIT_TAG v0.2.6
    PREFIX ${CMAKE_BINARY_DIR}/_deps/qoco-cpu
    GIT_SHALLOW TRUE
    UPDATE_DISCONNECTED TRUE
    CMAKE_ARGS ${QOCO_CMAKE_ARGS_CPU}
)

if(QOCO_HAVE_CUDA_BACKEND)
    ExternalProject_Add(qoco_project_cuda
        GIT_REPOSITORY https://github.com/qoco-org/qoco.git
        GIT_TAG v0.2.6
        PREFIX ${CMAKE_BINARY_DIR}/_deps/qoco-cuda
        GIT_SHALLOW TRUE
        UPDATE_DISCONNECTED TRUE
        CMAKE_ARGS ${QOCO_CMAKE_ARGS_CUDA}
    )
endif()

add_custom_target(qoco_project)
add_dependencies(qoco_project qoco_project_cpu)
if(TARGET qoco_project_cuda)
    add_dependencies(qoco_project qoco_project_cuda)
endif()

# Clarabel via FetchContent. This reuses Eigen3::Eigen from parent CMake and
# exposes libclarabel_c_shared for direct linking.
FetchContent_Declare(
    Clarabel_cpp
    GIT_REPOSITORY https://github.com/oxfordcontrol/Clarabel.cpp.git
    GIT_TAG main
    GIT_SHALLOW TRUE
    GIT_SUBMODULES_RECURSE TRUE
)
FetchContent_MakeAvailable(Clarabel_cpp)

# Add subdirs
add_subdirectory(src)
if(TARGET qoco_project)
    add_dependencies(cardillo qoco_project)
endif()
add_subdirectory(examples)
