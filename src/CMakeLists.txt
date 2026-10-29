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
# FetchContent_Declare(coal GIT_REPOSITORY https://github.com/coal-library/coal.git GIT_TAG v3.0.4)
# Coal fork 41pha1/coal (branch beam-support-set-hook, on top of devel incl. coal-library/coal#908):
# custom shapes (GEOM_CUSTOM) for the beam hull colliders.
FetchContent_Declare(coal GIT_REPOSITORY https://github.com/41pha1/coal.git GIT_TAG a90ca6e50bc6a038598148c69cbd84afd9c19df2)
FetchContent_Declare(qoco GIT_REPOSITORY https://github.com/qoco-org/qoco.git GIT_TAG v0.3.2)
# FetchContent_Declare(conicxx SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../ConicXX")
FetchContent_Declare(conicxx GIT_REPOSITORY https://github.com/ConicXX/ConicXX.git GIT_TAG main)
FetchContent_MakeAvailable(Eigen3)

add_library(eigen3 INTERFACE IMPORTED GLOBAL)
target_link_libraries(eigen3 INTERFACE eigen)
# Clarabel.cpp passes --release to cargo only if CMAKE_BUILD_TYPE matches "Release", so this
# project's default RelWithDebInfo (and MinSizeRel) would build Clarabel's Rust core in cargo's
# unoptimized debug profile (about 25-60x slower solves). Present the optimized build types as
# Release to its CMake; Debug keeps the debug profile.
set(CARDILLO_SAVED_BUILD_TYPE "${CMAKE_BUILD_TYPE}")
if(CMAKE_BUILD_TYPE MATCHES "^(RelWithDebInfo|MinSizeRel)$")
    set(CMAKE_BUILD_TYPE Release)
endif()
FetchContent_MakeAvailable(clarabel)
set(CMAKE_BUILD_TYPE "${CARDILLO_SAVED_BUILD_TYPE}")
FetchContent_MakeAvailable(entt qoco coal)
# Treat Coal's headers as system headers: Coal's own inline code (narrowphase.h) uses members it has
# deprecated itself, which would otherwise warn in every translation unit that includes Coal.
get_target_property(COAL_INCLUDE_DIRS coal INTERFACE_INCLUDE_DIRECTORIES)
set_target_properties(coal PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${COAL_INCLUDE_DIRS}")


FetchContent_MakeAvailable(conicxx)

# --- 2. Engine Targets ---
file(GLOB_RECURSE CARDILLO_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
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