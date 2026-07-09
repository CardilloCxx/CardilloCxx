cmake_minimum_required(VERSION 3.24)
project(cardillo LANGUAGES CXX)
set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "Install path prefix" FORCE)
include(FetchContent)
find_package(OpenMP QUIET)

# --- 1. Dependencies ---
set(BUILD_PYTHON_INTERFACE OFF CACHE BOOL "Disable Coal Python bindings" FORCE)
set(COAL_BUILD_PYTHON_INTERFACE OFF CACHE BOOL "Disable Coal Python bindings" FORCE)
set(COAL_BUILD_TESTING OFF CACHE BOOL "Disable Coal testing" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
# Eigen >= 5.0 gates its install(EXPORT Eigen3Targets ...) call behind this
# option, defaulted to PROJECT_IS_TOP_LEVEL (OFF for a FetchContent subproject).
# Without it, the "Eigen3Targets" export set is registered via install(TARGETS)
# but never completed via install(EXPORT), so coal's install(EXPORT coalTargets)
# fails dependency-completeness validation ("target eigen ... not in any export
# set"). Eigen 3.4 called install(EXPORT) unconditionally, so this never came up.
set(EIGEN_BUILD_CMAKE_PACKAGE ON CACHE BOOL "" FORCE)
set(CONICXX_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(CONICXX_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(CONICXX_BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)
set(CONICXX_INSTALL OFF CACHE BOOL "" FORCE)

FetchContent_Declare(eigen GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git GIT_TAG 5.0.1)
FetchContent_Declare(entt GIT_REPOSITORY https://github.com/skypjack/entt.git GIT_TAG v3.13.2)
FetchContent_Declare(clarabel GIT_REPOSITORY https://github.com/oxfordcontrol/Clarabel.cpp.git GIT_TAG v0.11.1)
FetchContent_Declare(coal GIT_REPOSITORY https://github.com/coal-library/coal.git GIT_TAG v3.0.0)
FetchContent_Declare(qoco GIT_REPOSITORY https://github.com/qoco-org/qoco.git GIT_TAG v0.3.2)
# Pinned past the initial "First working solver." integration point to
# cddaa5d ("Custom ldlt backend with dynamic regularization."), which adds
# two alternative KKT linear-solver backends (QDLDL, and a custom
# regularized-LDLT) on top of the "Speed fixes" (190eb67) and "Improved ill
# conditioned problems" (7dc4fab) commits already folded in at this point.
FetchContent_Declare(conicxx GIT_REPOSITORY https://github.com/ConicXX/ConicXX.git GIT_TAG cddaa5d5454fc0eb82e7b5851a09fd07cf6777f3)
FetchContent_MakeAvailable(eigen coal entt clarabel qoco)

# QOCO vendors its own copy of QDLDL under lib/qdldl (plain add_subdirectory,
# not FetchContent) at *global* CMake target scope: it unconditionally defines
# targets literally named `qdldlobject`/`qdldlstatic`. ConicXX independently
# wants a real, current qdldl (github.com/osqp/qdldl v0.1.9) as one of its
# selectable KKT linear-solver backends (see
# docs/chapters/solvers/interior_point.rst), fetched via its own
# `find_package(qdldl QUIET)` + FetchContent-fallback using those exact same
# unqualified target names.
#
# These two copies are NOT interchangeable: QOCO's vendored copy is an older
# qdldl whose QDLDL_factor() has a different (shorter) parameter list than
# the v0.1.9 API ConicXX's kkt/qdldl_ldlt.cpp and kkt/regularized_ldlt.cpp are
# written against -- aliasing one copy for both call sites was tried first
# and fails to compile with a hard "too few arguments to function
# QDLDL_factor" error. So both copies are genuinely needed simultaneously,
# which means their global CMake target names must not collide.
#
# FetchContent tracks population per *content name*, not per declared
# GIT_REPOSITORY/GIT_TAG, so if "qdldl" is already populated under that name
# by the time ConicXX's own (re-)declare + FetchContent_MakeAvailable(qdldl)
# call runs, that call is a no-op and ConicXX picks up whatever we already
# provided. So: fetch the real qdldl v0.1.9 ourselves first, with a
# PATCH_COMMAND that renames its `qdldlobject`/`qdldlstatic` targets to
# `conicxx_qdldlobject`/`conicxx_qdldlstatic` (see
# cmake/qdldl_rename_targets.sed) so they can't collide with QOCO's copy of
# the same names, then alias `qdldl::qdldl` to the renamed static target --
# exactly what ConicXX's own (unpatched) logic would have aliased, just
# under collision-safe names. (Flagged upstream to ConicXX too: an optional
# dependency that always claims global, unnamespaced CMake target names is
# inherently fragile against any sibling project vendoring the same
# library under an incompatible version -- QOCO isn't unique in doing this.)
set(QDLDL_BUILD_STATIC_LIB ON CACHE BOOL "" FORCE)
set(QDLDL_BUILD_SHARED_LIB OFF CACHE BOOL "" FORCE)
set(QDLDL_BUILD_DEMO_EXE OFF CACHE BOOL "" FORCE)
set(QDLDL_UNITTESTS OFF CACHE BOOL "" FORCE)
set(QDLDL_LONG OFF CACHE BOOL "" FORCE)  # QDLDL_int == int, matching conicxx::Index
FetchContent_Declare(qdldl
    GIT_REPOSITORY https://github.com/osqp/qdldl.git
    GIT_TAG v0.1.9
    GIT_SHALLOW TRUE
    PATCH_COMMAND sed -i -f "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/qdldl_rename_targets.sed" CMakeLists.txt
)
FetchContent_MakeAvailable(qdldl)
if(NOT TARGET qdldl::qdldl)
    add_library(qdldl::qdldl ALIAS conicxx_qdldlstatic)
endif()
FetchContent_MakeAvailable(conicxx)

# --- 2. Engine Targets ---
file(GLOB_RECURSE CARDILLO_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
add_library(cardillo SHARED ${CARDILLO_SOURCES})

# --- 3. Link Dependencies ---
target_link_libraries(cardillo PUBLIC
    Eigen3::Eigen
    EnTT::EnTT
    coal::coal
    qoco
    libclarabel_c_shared
    conicxx::conicxx
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(cardillo PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(cardillo PUBLIC CARDILLO_HAVE_OPENMP=1)
else()
    message(STATUS "OpenMP not found; CondensedSolver's jacobi/colored/chaotic sweep modes will run single-threaded.")
endif()

target_include_directories(cardillo PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

# --- 4. QOCO-Specific Integration ---
# Since QOCO requires specific path definitions, we derive them from the FetchContent locations
set(QOCO_CPU_LIB_DIR "${qoco_BINARY_DIR}/lib")
target_include_directories(cardillo PUBLIC 
    "${qoco_SOURCE_DIR}/include"
    "${qoco_SOURCE_DIR}/include/qdldl"
    "${qoco_SOURCE_DIR}/src/include"
)

get_target_property(QOCO_INCLUDES qoco INTERFACE_INCLUDE_DIRECTORIES)
if(NOT QOCO_INCLUDES)
    target_include_directories(qoco INTERFACE "${qoco_SOURCE_DIR}/include")
endif()

target_compile_definitions(cardillo PUBLIC 
    HAVE_QOCO=1
    QOCO_CPU_LIB_PATH="${QOCO_CPU_LIB_DIR}/libqoco.so"
    QOCO_CUDA_LIB_PATH="${QOCO_CPU_LIB_DIR}/libqoco.so" 
)

# --- 5. Clean Binary Properties ---
target_compile_definitions(cardillo PRIVATE BUILDING_CARDILLO_LIB=1)
# set_target_properties(cardillo PROPERTIES
#     CXX_VISIBILITY_PRESET hidden
#     VISIBILITY_INLINES_HIDDEN 1
#     POSITION_INDEPENDENT_CODE ON
# )