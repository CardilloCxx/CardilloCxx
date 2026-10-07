# FindMOSEK.cmake -- locate a user/system installation of the MOSEK Optimizer API (C, mosek.h).
#
# MOSEK is not redistributed with CardilloCxx; install it yourself (https://www.mosek.com/downloads/)
# and point this module at it with either
#   -DMOSEK_ROOT=<dir>         (CMake variable or environment variable), or
#   MOSEK_HOME=<dir>           (environment variable),
# where <dir> is the installation's platform directory (e.g. ~/mosek/11.2/tools/platform/linux64x86)
# or any directory above it (e.g. ~/mosek/11.2 or ~/mosek). Without hints, ~/mosek/<version> is searched.
#
# Defines the imported target MOSEK::MOSEK and MOSEK_FOUND, MOSEK_INCLUDE_DIR, MOSEK_LIBRARY,
# MOSEK_VERSION.

if(WIN32)
    set(_mosek_platform win64x86)
elseif(APPLE)
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")
        set(_mosek_platform osxaarch64)
    else()
        set(_mosek_platform osx64x86)
    endif()
else()
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64")
        set(_mosek_platform linuxaarch64)
    else()
        set(_mosek_platform linux64x86)
    endif()
endif()

set(_mosek_roots)
foreach(_r IN ITEMS "${MOSEK_ROOT}" "$ENV{MOSEK_ROOT}" "$ENV{MOSEK_HOME}")
    if(_r)
        list(APPEND _mosek_roots "${_r}")
    endif()
endforeach()
# Newest ~/mosek/<version> first.
file(GLOB _mosek_home_versions LIST_DIRECTORIES true "$ENV{HOME}/mosek/*")
list(SORT _mosek_home_versions COMPARE NATURAL ORDER DESCENDING)
list(APPEND _mosek_roots ${_mosek_home_versions})

set(_mosek_hints)
foreach(_r IN LISTS _mosek_roots)
    list(APPEND _mosek_hints "${_r}" "${_r}/tools/platform/${_mosek_platform}")
    file(GLOB _sub LIST_DIRECTORIES true "${_r}/*/tools/platform/${_mosek_platform}")
    list(SORT _sub COMPARE NATURAL ORDER DESCENDING)
    list(APPEND _mosek_hints ${_sub})
endforeach()

find_path(MOSEK_INCLUDE_DIR mosek.h HINTS ${_mosek_hints} PATH_SUFFIXES h include)

if(MOSEK_INCLUDE_DIR AND EXISTS "${MOSEK_INCLUDE_DIR}/mosek.h")
    file(STRINGS "${MOSEK_INCLUDE_DIR}/mosek.h" _mosek_ver REGEX "^#define MSK_VERSION_(MAJOR|MINOR|REVISION) +[0-9]+")
    string(REGEX REPLACE ".*MSK_VERSION_MAJOR +([0-9]+).*" "\\1" _mosek_major "${_mosek_ver}")
    string(REGEX REPLACE ".*MSK_VERSION_MINOR +([0-9]+).*" "\\1" _mosek_minor "${_mosek_ver}")
    string(REGEX REPLACE ".*MSK_VERSION_REVISION +([0-9]+).*" "\\1" _mosek_rev "${_mosek_ver}")
    set(MOSEK_VERSION "${_mosek_major}.${_mosek_minor}.${_mosek_rev}")
    get_filename_component(_mosek_platform_dir "${MOSEK_INCLUDE_DIR}" DIRECTORY)
endif()

find_library(MOSEK_LIBRARY
    NAMES mosek64 "mosek64_${_mosek_major}_${_mosek_minor}"
    HINTS "${_mosek_platform_dir}" ${_mosek_hints}
    PATH_SUFFIXES bin lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MOSEK
    REQUIRED_VARS MOSEK_LIBRARY MOSEK_INCLUDE_DIR
    VERSION_VAR MOSEK_VERSION)

if(MOSEK_FOUND AND NOT TARGET MOSEK::MOSEK)
    add_library(MOSEK::MOSEK UNKNOWN IMPORTED)
    set_target_properties(MOSEK::MOSEK PROPERTIES
        IMPORTED_LOCATION "${MOSEK_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${MOSEK_INCLUDE_DIR}")
endif()

mark_as_advanced(MOSEK_INCLUDE_DIR MOSEK_LIBRARY)
