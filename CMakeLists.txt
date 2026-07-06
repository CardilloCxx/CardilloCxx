# 3.24 is required for FetchContent's OVERRIDE_FIND_PACKAGE/FetchContent_SetPopulated,
# used in src/CMakeLists.txt to deduplicate QOCO's vendored qdldl copy against
# ConicXX's own qdldl dependency (see the comment there for the full story).
cmake_minimum_required(VERSION 3.24)
project(cardillo_project LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# CMake's own default build type is empty (no -O/-DNDEBUG at all) unless explicitly requested;
# for a numerics-heavy library that default silently produces an unoptimized build. Default to
# RelWithDebInfo (-O2 -g -DNDEBUG) so a plain `cmake ..` is fast, while still allowing
# `-DCMAKE_BUILD_TYPE=...` to override it as usual. Skipped for multi-config generators
# (Visual Studio/Xcode), which select the configuration at build time instead.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
endif()

add_subdirectory(src)
add_subdirectory(examples)