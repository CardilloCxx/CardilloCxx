# Finds the C API of NVIDIA cuOpt (library cuopt_mathopt, header cuopt/mathematical_optimization/cuopt_c.h).
# cuOpt's own CMake package needs CMake >= 4.0 and is therefore not used.
#
# Searched below CUOPT_ROOT (CMake or environment variable), VIRTUAL_ENV and CONDA_PREFIX. Each may be a
# conda prefix, a Python virtual environment (pip package libcuopt-mathopt-cu12/-cu13) or its
# site-packages directory. pip splits the headers: cuopt_c.h lies in libcuopt_mathopt, the headers it
# includes (constants.h, export.hpp, version_config.hpp) in libcuopt_client.
#
# Defines CUOPT_FOUND, CUOPT_VERSION, CUOPT_INCLUDE_DIRS, CUOPT_LIBRARY and the imported target CUOPT::CUOPT.

set(_cuopt_roots)
foreach(_r IN ITEMS "${CUOPT_ROOT}" "$ENV{CUOPT_ROOT}" "$ENV{VIRTUAL_ENV}" "$ENV{CONDA_PREFIX}")
    if(_r)
        list(APPEND _cuopt_roots "${_r}")
    endif()
endforeach()

set(_cuopt_hints)
foreach(_r IN LISTS _cuopt_roots)
    list(APPEND _cuopt_hints "${_r}" "${_r}/libcuopt_mathopt" "${_r}/libcuopt_client")
    file(GLOB _site LIST_DIRECTORIES true "${_r}/lib/python3*/site-packages")
    foreach(_s IN LISTS _site)
        list(APPEND _cuopt_hints "${_s}/libcuopt_mathopt" "${_s}/libcuopt_client")
    endforeach()
endforeach()

find_path(CUOPT_INCLUDE_DIR cuopt/mathematical_optimization/cuopt_c.h HINTS ${_cuopt_hints} PATH_SUFFIXES include)
find_path(CUOPT_CLIENT_INCLUDE_DIR cuopt/version_config.hpp HINTS ${_cuopt_hints} PATH_SUFFIXES include)
find_library(CUOPT_LIBRARY NAMES cuopt_mathopt HINTS ${_cuopt_hints} PATH_SUFFIXES lib64 lib)

if(CUOPT_CLIENT_INCLUDE_DIR AND EXISTS "${CUOPT_CLIENT_INCLUDE_DIR}/cuopt/version_config.hpp")
    file(STRINGS "${CUOPT_CLIENT_INCLUDE_DIR}/cuopt/version_config.hpp" _cuopt_ver REGEX "^#define CUOPT_VERSION_(MAJOR|MINOR|PATCH) +[0-9]+")
    string(REGEX REPLACE ".*CUOPT_VERSION_MAJOR +([0-9]+).*" "\\1" _cuopt_major "${_cuopt_ver}")
    string(REGEX REPLACE ".*CUOPT_VERSION_MINOR +([0-9]+).*" "\\1" _cuopt_minor "${_cuopt_ver}")
    string(REGEX REPLACE ".*CUOPT_VERSION_PATCH +([0-9]+).*" "\\1" _cuopt_patch "${_cuopt_ver}")
    set(CUOPT_VERSION "${_cuopt_major}.${_cuopt_minor}.${_cuopt_patch}")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CUOPT
    REQUIRED_VARS CUOPT_LIBRARY CUOPT_INCLUDE_DIR CUOPT_CLIENT_INCLUDE_DIR
    VERSION_VAR CUOPT_VERSION)

if(CUOPT_FOUND)
    set(CUOPT_INCLUDE_DIRS "${CUOPT_INCLUDE_DIR}" "${CUOPT_CLIENT_INCLUDE_DIR}")
    list(REMOVE_DUPLICATES CUOPT_INCLUDE_DIRS)
    if(NOT TARGET CUOPT::CUOPT)
        add_library(CUOPT::CUOPT UNKNOWN IMPORTED)
        set_target_properties(CUOPT::CUOPT PROPERTIES
            IMPORTED_LOCATION "${CUOPT_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${CUOPT_INCLUDE_DIRS}")
    endif()
endif()

mark_as_advanced(CUOPT_INCLUDE_DIR CUOPT_CLIENT_INCLUDE_DIR CUOPT_LIBRARY)
