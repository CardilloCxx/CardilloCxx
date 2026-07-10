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

# ThreadSanitizer instrumentation, for validating CondensedSolver's parallel (jacobi/colored/
# chaotic) sweep modes. TSan needs -O1 (not -O0, and not RelWithDebInfo's -O2 -DNDEBUG) for a
# reasonable instrumentation/perf tradeoff with usable stack traces, so configure with
# -DCARDILLO_ENABLE_TSAN=ON -DCMAKE_BUILD_TYPE=Debug (or any non-RelWithDebInfo type) -- otherwise
# RelWithDebInfo's own -O2 -DNDEBUG fights these flags. Applied at this top-level scope, before
# add_subdirectory, so it also instruments third-party C/C++ deps fetched under src/ (they're
# add_subdirectory'd from there); Clarabel's Rust core is a separate cargo build and is NOT
# instrumented -- TSan won't observe races inside it, which is fine since CondensedSolver never
# calls into it.
option(CARDILLO_ENABLE_TSAN "Build with ThreadSanitizer instrumentation" OFF)
if(CARDILLO_ENABLE_TSAN)
    add_compile_options(-fsanitize=thread -g -O1)
    add_link_options(-fsanitize=thread)
endif()

add_subdirectory(src)
add_subdirectory(examples)
# add_subdirectory(frontend)