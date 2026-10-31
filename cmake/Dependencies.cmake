# Required dependencies. Eigen and EnTT are taken from an installation (or a parent project) if
# find_package() finds them, otherwise they are fetched. Coal is the fork 41pha1/coal (see below).
include(FetchContent)

# Settings for the fetched dependencies are normal variables of this directory scope, so they do not
# leak into the cache of a parent project; CMP0077 makes option() in the dependencies honor them.
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
set(BUILD_TESTING OFF)
set(BUILD_EXAMPLES OFF)
set(BUILD_PYTHON_INTERFACE OFF)
set(COAL_BUILD_PYTHON_INTERFACE OFF)
set(COAL_BUILD_TESTING OFF)
set(EIGEN_BUILD_DOC OFF)
set(EIGEN_BUILD_TESTING OFF)
# Eigen >= 5.0 installs its CMake package only when this is on (default: top-level project only);
# Coal's and cardillo's exported targets need it.
set(EIGEN_BUILD_CMAKE_PACKAGE ON)

find_package(OpenMP QUIET)

# --- VTK (VTK export, private to the library) ---
# Prefer a system VTK. Debian/Ubuntu ParaView packages ship their own VTK under a non-standard
# location (and conflict with libvtk9-dev), so probe that as well before building VTK from source.
set(CARDILLO_VTK_COMPONENTS CommonCore CommonDataModel IOXML)
find_package(VTK QUIET COMPONENTS ${CARDILLO_VTK_COMPONENTS})
if(NOT VTK_FOUND)
    file(GLOB CARDILLO_VTK_PARAVIEW_HINTS "/usr/lib/*/cmake/paraview-*/vtk" "/usr/local/lib/*/cmake/paraview-*/vtk")
    if(CARDILLO_VTK_PARAVIEW_HINTS)
        find_package(VTK QUIET COMPONENTS ${CARDILLO_VTK_COMPONENTS} PATHS ${CARDILLO_VTK_PARAVIEW_HINTS} NO_DEFAULT_PATH)
    endif()
endif()
if(VTK_FOUND)
    message(STATUS "Using system VTK ${VTK_VERSION} (${VTK_DIR})")
else()
    message(STATUS "No system VTK found; building VTK from source (slow)")
    set(VTK_GROUP_ENABLE_StandAlone "DONT_WANT")
    set(VTK_MODULE_ENABLE_VTK_CommonCore "YES")
    set(VTK_MODULE_ENABLE_VTK_CommonDataModel "YES")
    set(VTK_MODULE_ENABLE_VTK_IOXML "YES")
    set(VTK_BUILD_TESTING OFF)
    set(VTK_BUILD_EXAMPLES OFF)
    set(VTK_BUILD_DOCUMENTATION OFF)
    set(VTK_WRAP_PYTHON OFF)
    set(VTK_WRAP_JAVA OFF)
    FetchContent_Declare(vtk GIT_REPOSITORY https://gitlab.kitware.com/vtk/vtk.git GIT_TAG v9.3.1)
    FetchContent_MakeAvailable(vtk)
endif()

# --- Eigen, EnTT ---
# Eigen's cmake_dependent_option(BUILD_TESTING ...) force-sets the cache entry to OFF when Eigen is not
# the top-level project; a parent project's value is restored after Coal (the last dependency).
if(DEFINED CACHE{BUILD_TESTING})
    set(CARDILLO_PARENT_BUILD_TESTING "$CACHE{BUILD_TESTING}")
endif()
FetchContent_Declare(Eigen3 GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git GIT_TAG 5.0.1 FIND_PACKAGE_ARGS 5.0 NAMES Eigen3)
FetchContent_Declare(EnTT GIT_REPOSITORY https://github.com/skypjack/entt.git GIT_TAG v3.13.2 FIND_PACKAGE_ARGS 3.13 NAMES EnTT)
FetchContent_MakeAvailable(Eigen3 EnTT)

# --- Coal ---
# The fork 41pha1/coal (branch beam-support-set-hook, on top of devel incl. coal-library/coal#908) adds
# custom shapes (GEOM_CUSTOM) for the beam hull colliders; released Coal versions do not build with
# CardilloCxx. See .claude/docs/collision.md.
if(CARDILLO_USE_SYSTEM_COAL)
    find_package(coal REQUIRED)
else()
    FetchContent_Declare(coal GIT_REPOSITORY https://github.com/41pha1/coal.git GIT_TAG a90ca6e50bc6a038598148c69cbd84afd9c19df2)
    FetchContent_MakeAvailable(coal)
    # Coal's own inline code uses members it has deprecated itself; treat its headers as system headers.
    get_target_property(_coal_include_dirs coal INTERFACE_INCLUDE_DIRECTORIES)
    set_target_properties(coal PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${_coal_include_dirs}")
    # Warnings from compiling Coal itself are not ours to fix.
    target_compile_options(coal PRIVATE -w)
endif()

if(DEFINED CARDILLO_PARENT_BUILD_TESTING)
    set_property(CACHE BUILD_TESTING PROPERTY VALUE "${CARDILLO_PARENT_BUILD_TESTING}")
endif()
