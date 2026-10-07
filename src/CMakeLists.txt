cmake_minimum_required(VERSION 3.24)
project(cardillo LANGUAGES CXX)
set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "Install path prefix" FORCE)
include(FetchContent)
include(ExternalProject)
find_package(OpenMP QUIET)

# --- VTK (used by the VTK export writer) ---
# Prefer a system-installed VTK (fast, no compile) over building it from source.
# 1) Standard find_package: covers a normal libvtk9-dev-style install.
# 2) Debian/Ubuntu's ParaView packages bundle their own VTK under a non-standard
#    CMake config location (/usr/lib/*/cmake/paraview-*/vtk instead of the usual
#    vtk-9.x path) -- installing the standalone libvtk9-dev package conflicts with
#    an existing ParaView install on such systems (it depends on python3-vtk9,
#    which Conflicts: python3-paraview), so probe that location as a fallback.
# 3) Only build VTK from source via FetchContent if neither is found.
set(CARDILLO_VTK_COMPONENTS CommonCore CommonDataModel IOXML)
find_package(VTK QUIET COMPONENTS ${CARDILLO_VTK_COMPONENTS})
if(NOT VTK_FOUND)
    file(GLOB CARDILLO_VTK_PARAVIEW_HINTS
        "/usr/lib/*/cmake/paraview-*/vtk"
        "/usr/local/lib/*/cmake/paraview-*/vtk"
    )
    if(CARDILLO_VTK_PARAVIEW_HINTS)
        find_package(VTK QUIET COMPONENTS ${CARDILLO_VTK_COMPONENTS} PATHS ${CARDILLO_VTK_PARAVIEW_HINTS} NO_DEFAULT_PATH)
    endif()
endif()

if(VTK_FOUND)
    message(STATUS "Using system VTK ${VTK_VERSION} (${VTK_DIR}) for VTK export")
else()
    message(STATUS "No system VTK found; building VTK from source via FetchContent (this is slow)")
    set(VTK_GROUP_ENABLE_StandAlone "DONT_WANT" CACHE STRING "" FORCE)
    set(VTK_MODULE_ENABLE_VTK_CommonCore "YES" CACHE STRING "" FORCE)
    set(VTK_MODULE_ENABLE_VTK_CommonDataModel "YES" CACHE STRING "" FORCE)
    set(VTK_MODULE_ENABLE_VTK_IOXML "YES" CACHE STRING "" FORCE)
    set(VTK_BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(VTK_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(VTK_BUILD_DOCUMENTATION OFF CACHE BOOL "" FORCE)
    set(VTK_WRAP_PYTHON OFF CACHE BOOL "" FORCE)
    set(VTK_WRAP_JAVA OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(vtk GIT_REPOSITORY https://gitlab.kitware.com/vtk/vtk.git GIT_TAG v9.3.1)
    FetchContent_MakeAvailable(vtk)
endif()

# --- 0. Set paths ---
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/config/path.hpp.in
    ${CMAKE_CURRENT_SOURCE_DIR}/config/path.hpp
    @ONLY
)

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

FetchContent_Declare(Eigen3 GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git GIT_TAG 5.0.1 OVERRIDE_FIND_PACKAGE)
FetchContent_Declare(entt GIT_REPOSITORY https://github.com/skypjack/entt.git GIT_TAG v3.13.2)
FetchContent_Declare(clarabel GIT_REPOSITORY https://github.com/oxfordcontrol/Clarabel.cpp.git GIT_TAG v0.11.1)
FetchContent_Declare(coal GIT_REPOSITORY https://github.com/coal-library/coal.git GIT_TAG v3.0.4)
FetchContent_Declare(qoco GIT_REPOSITORY https://github.com/qoco-org/qoco.git GIT_TAG v0.3.2)
FetchContent_Declare(conicxx GIT_REPOSITORY https://github.com/ConicXX/ConicXX.git GIT_TAG cddaa5d5454fc0eb82e7b5851a09fd07cf6777f3)
FetchContent_MakeAvailable(Eigen3)

add_library(eigen3 INTERFACE IMPORTED GLOBAL)
target_link_libraries(eigen3 INTERFACE eigen)
# Clarabel.cpp passes --release to cargo only if CMAKE_BUILD_TYPE matches "Release", so this
# project's default RelWithDebInfo (and MinSizeRel) silently built Clarabel's Rust core in cargo's
# unoptimized debug profile (~25-60x slower solves). Present optimized C/C++ build types as Release
# to its CMake; Debug keeps the debug profile.
set(CARDILLO_SAVED_BUILD_TYPE "${CMAKE_BUILD_TYPE}")
if(CMAKE_BUILD_TYPE MATCHES "^(RelWithDebInfo|MinSizeRel)$")
    set(CMAKE_BUILD_TYPE Release)
endif()
FetchContent_MakeAvailable(clarabel)
set(CMAKE_BUILD_TYPE "${CARDILLO_SAVED_BUILD_TYPE}")
FetchContent_MakeAvailable(entt qoco coal)


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
    PATCH_COMMAND sed -i.bak -f "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/qdldl_rename_targets.sed" CMakeLists.txt
)
FetchContent_MakeAvailable(qdldl)
if(NOT TARGET qdldl::qdldl)
    add_library(qdldl::qdldl ALIAS conicxx_qdldlstatic)
endif()
FetchContent_MakeAvailable(conicxx)

# --- 2. Engine Targets ---
file(GLOB_RECURSE CARDILLO_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
# MOSEK backend sources are only compiled when the optional dependency is enabled, so no MOSEK
# header or symbol is referenced otherwise.
if(NOT CARDILLO_WITH_MOSEK)
    list(FILTER CARDILLO_SOURCES EXCLUDE REGEX "/physics/solver/mosek_[^/]*\\.cpp$")
endif()
if(NOT CARDILLO_WITH_SCS)
    list(FILTER CARDILLO_SOURCES EXCLUDE REGEX "/physics/solver/scs/[^/]*\\.cpp$")
endif()
add_library(cardillo SHARED ${CARDILLO_SOURCES})


# Python version of Coal:
# find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
# include_directories( ${PYTHON3_INCLUDE_DIRS})
# set(COAL_INSTALL_DIR ${CMAKE_BINARY_DIR}/external/coal)
# ExternalProject_Add(coal_project
#     PREFIX ${CMAKE_BINARY_DIR}/_deps/coal
#     DOWNLOAD_COMMAND ""
#     CONFIGURE_COMMAND ""
#     BUILD_COMMAND ${Python3_EXECUTABLE} -m pip install coal -t ${COAL_INSTALL_DIR}
#     INSTALL_COMMAND ""
# )
# 
# add_dependencies(${PROJECT_NAME} coal_project)
# target_include_directories(${PROJECT_NAME} PUBLIC
#     "${COAL_INSTALL_DIR}/cmeel.prefix/include/"
# )
# 
# target_link_directories(${PROJECT_NAME} PUBLIC
#     "${COAL_INSTALL_DIR}/cmeel.prefix/lib/"
# )
# target_link_libraries(${PROJECT_NAME} INTERFACE coal)


# --- 3. Link Dependencies ---
target_link_libraries(cardillo PUBLIC
    Eigen3::Eigen
    EnTT::EnTT
    coal::coal
    qoco
    libclarabel_c_shared
    conicxx::conicxx
    VTK::CommonCore
    VTK::CommonDataModel
    VTK::IOXML
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(cardillo PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(cardillo PUBLIC CARDILLO_HAVE_OPENMP=1)
else()
    message(STATUS "OpenMP not found; CondensedSolver's jacobi/colored/chaotic sweep modes will run single-threaded.")
endif()

target_include_directories(cardillo PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

# --- MOSEK (optional) ---
if(CARDILLO_WITH_MOSEK)
    find_package(MOSEK REQUIRED)
    message(STATUS "Using MOSEK ${MOSEK_VERSION} (${MOSEK_LIBRARY}) for the mosek solver backend")
    target_link_libraries(cardillo PUBLIC MOSEK::MOSEK)
    target_compile_definitions(cardillo PUBLIC CARDILLO_HAVE_MOSEK=1)
endif()

# --- SCS (optional) ---
# Only the direct backend (scsdir: AMD ordering + QDLDL LDL' of the quasi-definite KKT matrix) is
# used. SCS vendors its own AMD and QDLDL under their original, unprefixed symbol names, and both
# already exist in this process: libcardillo exports ConicXX's QDLDL (QDLDL_factor, ...) and libqoco
# exports QOCO's AMD (amd_order, ...). scsdir is therefore built as its own shared library that
# exports nothing but the public scs_* API (cmake/scs_exports.map): its internal QDLDL/AMD calls bind
# to its own copies, and it cannot interpose on ConicXX's or QOCO's. scs_int = int and
# scs_float = double (DLONG/SFLOAT off) are additionally checked by static_asserts in the adapter.
if(CARDILLO_WITH_SCS)
    set(BUILD_SHARED_LIBS ON)
    # SCS compiles Anderson acceleration only with USE_LAPACK (src/aa.c); without it AA is silently
    # disabled, which costs ADMM a large factor in iterations. Any BLAS/LAPACK will do (AA only needs
    # small dense least-squares solves).
    find_package(LAPACK QUIET)
    if(LAPACK_FOUND)
        set(USE_LAPACK ON)
    else()
        set(USE_LAPACK OFF)
        message(WARNING "CARDILLO_WITH_SCS: no BLAS/LAPACK found; SCS is built without Anderson acceleration (expect many more iterations)")
    endif()
    set(DLONG OFF)
    set(SFLOAT OFF)
    # SCS's cmake/AddUninstallTarget.cmake defines a global `uninstall` target, which coal already
    # defines; its own include guard skips it.
    set(__ADD_UNINSTALL_TARGET_INCLUDED TRUE)
    # SCS links OpenMP::OpenMP_C (and so compiles its `omp parallel for` loops over the cones and in
    # a sparse mat-vec) whenever OpenMP_C_FOUND is set, independent of its own USE_OPENMP, and our
    # find_package(OpenMP) above sets it. Fork-joining a full thread team per ADMM iteration over
    # 3-dimensional cones made small scenes ~7x slower in practice, so SCS is built serial unless
    # CARDILLO_SCS_OPENMP is on.
    set(CARDILLO_SAVED_OpenMP_C_FOUND "${OpenMP_C_FOUND}")
    if(NOT CARDILLO_SCS_OPENMP)
        set(OpenMP_C_FOUND FALSE)
    endif()
    FetchContent_Declare(scs GIT_REPOSITORY https://github.com/cvxgrp/scs.git GIT_TAG 3.3.1 GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(scs)
    unset(BUILD_SHARED_LIBS)
    unset(USE_LAPACK)
    unset(DLONG)
    unset(SFLOAT)
    unset(__ADD_UNINSTALL_TARGET_INCLUDED)
    set(OpenMP_C_FOUND "${CARDILLO_SAVED_OpenMP_C_FOUND}")
    # SCS always defines its indirect (CG) and, with LAPACK, dense libraries; only scsdir is used.
    foreach(_scs_unused scsindir scsdense)
        if(TARGET ${_scs_unused})
            set_target_properties(${_scs_unused} PROPERTIES EXCLUDE_FROM_ALL ON)
        endif()
    endforeach()
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        target_link_options(scsdir PRIVATE "LINKER:--version-script=${CMAKE_CURRENT_SOURCE_DIR}/../cmake/scs_exports.map")
    endif()
    # (macOS: two-level namespaces already bind scsdir's internal references to its own copies.)
    target_link_libraries(cardillo PUBLIC scs::scsdir)
    target_compile_definitions(cardillo PUBLIC CARDILLO_HAVE_SCS=1)
endif()

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
    # The shared library's actual build location (it is not under ${QOCO_CPU_LIB_DIR}).
    QOCO_CPU_LIB_PATH="$<TARGET_FILE:qoco>"
    QOCO_CUDA_LIB_PATH="$<TARGET_FILE:qoco>" 
)

# --- 5. Clean Binary Properties ---
target_compile_definitions(cardillo PRIVATE BUILDING_CARDILLO_LIB=1)
# set_target_properties(cardillo PROPERTIES
#     CXX_VISIBILITY_PRESET hidden
#     VISIBILITY_INLINES_HIDDEN 1
#     POSITION_INDEPENDENT_CODE ON
# )