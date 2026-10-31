# Optional conic solver backends, each guarded by its CARDILLO_WITH_<NAME> option. The sources of a
# backend live in src/cardillo/physics/solver/conic/<name>/ and are added to the library in
# src/CMakeLists.txt.
include(FetchContent)

if(CARDILLO_WITH_CLARABEL)
    FetchContent_Declare(clarabel GIT_REPOSITORY https://github.com/oxfordcontrol/Clarabel.cpp.git GIT_TAG v0.11.1)
    # Clarabel.cpp passes --release to cargo only for CMAKE_BUILD_TYPE=Release; present the other
    # optimized build types as Release so its Rust core is not built in the (25-60x slower) debug profile.
    set(CARDILLO_SAVED_BUILD_TYPE "${CMAKE_BUILD_TYPE}")
    if(CMAKE_BUILD_TYPE MATCHES "^(RelWithDebInfo|MinSizeRel)$")
        set(CMAKE_BUILD_TYPE Release)
    endif()
    FetchContent_MakeAvailable(clarabel)
    set(CMAKE_BUILD_TYPE "${CARDILLO_SAVED_BUILD_TYPE}")
endif()

if(CARDILLO_WITH_QOCO)
    FetchContent_Declare(qoco GIT_REPOSITORY https://github.com/qoco-org/qoco.git GIT_TAG v0.3.2)
    FetchContent_MakeAvailable(qoco)
endif()

if(CARDILLO_WITH_CONICXX)
    set(CONICXX_BUILD_TESTS OFF)
    set(CONICXX_BUILD_EXAMPLES OFF)
    set(CONICXX_BUILD_BENCHMARKS OFF)
    set(CONICXX_INSTALL OFF)  # static library, linked into libcardillo
    # A local checkout can be used with -DFETCHCONTENT_SOURCE_DIR_CONICXX=<path>.
    FetchContent_Declare(conicxx GIT_REPOSITORY https://github.com/ConicXX/ConicXX.git GIT_TAG main)
    FetchContent_MakeAvailable(conicxx)
endif()

# SCS: only its direct linear-system backend (scsdir: AMD + QDLDL) is used. SCS vendors AMD and QDLDL
# under their original symbol names, which may also exist in this process (ConicXX's QDLDL in
# libcardillo, QOCO's AMD in libqoco). scsdir is therefore a shared library that exports only the
# public scs_* API (cmake/scs_exports.map), so its internal calls bind to its own copies.
if(CARDILLO_WITH_SCS)
    # SCS compiles Anderson acceleration only with LAPACK; without it ADMM needs many more iterations.
    find_package(LAPACK QUIET)
    if(LAPACK_FOUND)
        set(USE_LAPACK ON)
    else()
        set(USE_LAPACK OFF)
        message(WARNING "CARDILLO_WITH_SCS: no BLAS/LAPACK found; SCS is built without Anderson acceleration")
    endif()
    set(BUILD_SHARED_LIBS ON)
    set(DLONG OFF)
    set(SFLOAT OFF)
    # SCS defines a global `uninstall` target, which Coal already defines; SCS's include guard skips it.
    set(__ADD_UNINSTALL_TARGET_INCLUDED TRUE)
    # SCS compiles its OpenMP loops whenever OpenMP_C_FOUND is set, independent of its own USE_OPENMP;
    # CARDILLO_SCS_OPENMP=OFF builds it serial (fork-joining per ADMM iteration can cost more than it gains).
    set(CARDILLO_SAVED_OpenMP_C_FOUND "${OpenMP_C_FOUND}")
    if(NOT CARDILLO_SCS_OPENMP)
        set(OpenMP_C_FOUND FALSE)
    endif()
    FetchContent_Declare(scs GIT_REPOSITORY https://github.com/cvxgrp/scs.git GIT_TAG 3.3.1 GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(scs)
    set(OpenMP_C_FOUND "${CARDILLO_SAVED_OpenMP_C_FOUND}")
    unset(BUILD_SHARED_LIBS)
    unset(USE_LAPACK)
    unset(DLONG)
    unset(SFLOAT)
    unset(__ADD_UNINSTALL_TARGET_INCLUDED)
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        target_link_options(scsdir PRIVATE "LINKER:--version-script=${PROJECT_SOURCE_DIR}/cmake/scs_exports.map")
    endif()
endif()

if(CARDILLO_WITH_MOSEK)
    find_package(MOSEK REQUIRED)
    message(STATUS "Using MOSEK ${MOSEK_VERSION} (${MOSEK_LIBRARY})")
endif()
