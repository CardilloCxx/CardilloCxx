cmake_minimum_required(VERSION 3.24)
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
project(FilamentExample LANGUAGES C CXX) # Added C here because Filament tests both compilers

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)

option(FILAMENT_ENABLE_EXPERIMENTAL_GCC_SUPPORT "Enable experimental GCC support for Filament" ON)
option(FILAMENT_BUILD_TESTING "Build Filament tests" OFF)
option(FILAMENT_SKIP_SAMPLES "Skip Filament samples" ON)
option(FILAMENT_SKIP_SDL2 "Skip SDL2-related Filament dependencies" ON)

FetchContent_Declare(
    filament
    GIT_REPOSITORY https://github.com/google/filament.git
    GIT_TAG main
    GIT_SHALLOW TRUE
)

FetchContent_Declare(
    glfw
    GIT_REPOSITORY https://github.com/glfw/glfw.git
    GIT_TAG 3.4
    GIT_SHALLOW TRUE
)

set(FILAMENT_ENABLE_EXPERIMENTAL_GCC_SUPPORT ${FILAMENT_ENABLE_EXPERIMENTAL_GCC_SUPPORT} CACHE BOOL "" FORCE)
set(FILAMENT_BUILD_TESTING ${FILAMENT_BUILD_TESTING} CACHE BOOL "" FORCE)
set(FILAMENT_SKIP_SAMPLES ${FILAMENT_SKIP_SAMPLES} CACHE BOOL "" FORCE)
set(FILAMENT_SKIP_SDL2 ${FILAMENT_SKIP_SDL2} CACHE BOOL "" FORCE)
set(FILAMENT_SUPPORTS_XCB ON CACHE BOOL "" FORCE)
set(FILAMENT_SUPPORTS_XLIB ON CACHE BOOL "" FORCE)
set(FILAMENT_SUPPORTS_WAYLAND OFF CACHE BOOL "" FORCE)
set(FILAMENT_SUPPORTS_EGL_ON_LINUX OFF CACHE BOOL "" FORCE)
set(FILAMENT_SUPPORTS_OSMESA OFF CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_WAYLAND OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_X11 ON CACHE BOOL "" FORCE)
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)

if(MSVC)
    add_compile_options(/wd4996) # The MSVC way to ignore deprecated warnings globally
else()
    add_compile_options(-Wno-error=deprecated-declarations)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-error=deprecated-declarations")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-error=deprecated-declarations")

    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-error=deprecated-declarations -stdlib=libc++")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-error=deprecated-declarations")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libc++ -lc++abi")
    endif()

    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        add_compile_options(-fpermissive)
    endif()
endif()
FetchContent_MakeAvailable(filament glfw)

find_package(OpenGL REQUIRED)

if(NOT WIN32)
    find_package(X11 REQUIRED)
endif()
set(ALL_FILAMATS "")

function(compile_filament_material MAT_NAME)
    set(MAT_SRC "${CMAKE_CURRENT_SOURCE_DIR}/assets/materials/${MAT_NAME}.mat")
    set(MAT_OUT "${CMAKE_CURRENT_BINARY_DIR}/generated/${MAT_NAME}.filamat")

    add_custom_command(
        OUTPUT "${MAT_OUT}"
        COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_BINARY_DIR}/generated"
        COMMAND matc -a vulkan -o "${MAT_OUT}" "${MAT_SRC}"
        DEPENDS "${MAT_SRC}"
        COMMENT "Compiling material: ${MAT_NAME}"
    )
    
    # Safely append to the global list using the current parent state
    set(TEMP_LIST ${ALL_FILAMATS})
    list(APPEND TEMP_LIST "${MAT_OUT}")
    set(ALL_FILAMATS "${TEMP_LIST}" PARENT_SCOPE)
endfunction()

compile_filament_material(pbr)
compile_filament_material(glass)
compile_filament_material(picking)

add_custom_target(filament_materials DEPENDS ${ALL_FILAMATS})

# ==============================================================================
# Texture Compilation
# ==============================================================================

# for f in albedo/*.jpg normal/*.jpg roughness/*.jpg  albedo/*.exr normal/*.exr roughness/*.exr; do
#   convert "$f" -resize "2048x2048>" "${f%.jpg}.png"
#   rm "$f"
# done

function(get_mipgen_flags REL_PATH OUT_FLAGS)
    get_filename_component(REL_DIR "${REL_PATH}" DIRECTORY)
    string(TOLOWER "${REL_DIR}" DIR_LOWER)

    set(FLAGS "")

    if(DIR_LOWER MATCHES "(^|/)normal($|/)")
        set(FLAGS --strip-alpha --kernel=NORMALS --linear)
    elseif(DIR_LOWER MATCHES "(^|/)(roughness|metallic|metalness|ao|occlusion)($|/)")
        set(FLAGS --grayscale --linear)
    elseif(DIR_LOWER MATCHES "(^|/)(albedo|diffuse|basecolor)($|/)")
        set(FLAGS "")
    else()
        message(STATUS "Texture '${REL_PATH}' is in an uncategorized folder. Using default sRGB mipgen settings.")
        set(FLAGS "")
    endif()

    set(${OUT_FLAGS} "${FLAGS}" PARENT_SCOPE)
endfunction()

file(GLOB_RECURSE TEXTURE_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/assets/textures/*.png"
    "${CMAKE_CURRENT_SOURCE_DIR}/assets/textures/*.jpg"
    "${CMAKE_CURRENT_SOURCE_DIR}/assets/textures/*.jpeg"
    "${CMAKE_CURRENT_SOURCE_DIR}/assets/textures/*.tga"
)

set(ALL_KTX_FILES "")
set(HEADER_ENTRIES "")

# Ensure the root generated textures directory exists during CMake configure phase
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/generated/textures")
foreach(TEX_PATH ${TEXTURE_SOURCES})
    file(RELATIVE_PATH REL_PATH "${CMAKE_CURRENT_SOURCE_DIR}/assets/textures" "${TEX_PATH}")
    get_filename_component(REL_DIR "${REL_PATH}" DIRECTORY)
    get_filename_component(FILE_NAME_WE "${REL_PATH}" NAME_WE)

    if(REL_DIR)
        set(REL_KTX_PATH "${REL_DIR}/${FILE_NAME_WE}.ktx")
    else()
        set(REL_KTX_PATH "${FILE_NAME_WE}.ktx")
    endif()

    set(KTX_OUT "${CMAKE_CURRENT_BINARY_DIR}/generated/textures/${REL_KTX_PATH}")
    get_filename_component(KTX_OUT_DIR "${KTX_OUT}" DIRECTORY)

    file(MAKE_DIRECTORY "${KTX_OUT_DIR}")

    get_mipgen_flags("${REL_PATH}" MIPGEN_FLAGS)

    add_custom_command(
        OUTPUT "${KTX_OUT}"
        COMMAND ${CMAKE_COMMAND} -E make_directory "${KTX_OUT_DIR}"
        COMMAND mipgen ${MIPGEN_FLAGS} "${TEX_PATH}" "${KTX_OUT}"
        DEPENDS "${TEX_PATH}" mipgen
        COMMENT "Compiling KTX [${REL_DIR}]: ${FILE_NAME_WE}"
        VERBATIM
    )

    list(APPEND ALL_KTX_FILES "${KTX_OUT}")

    string(REGEX REPLACE "[^a-zA-Z0-9_]" "_" CONST_NAME "${REL_KTX_PATH}")
    string(TOUPPER "${CONST_NAME}" CONST_NAME)
    
    string(APPEND HEADER_ENTRIES "    inline constexpr std::string_view ${CONST_NAME} = \"${REL_KTX_PATH}\";\n")
endforeach()

set(HEADER_CONTENT "#pragma once\n\n#include <string_view>\n\nnamespace Textures {\n    inline constexpr std::string_view BASE_DIR = \"generated/textures\";\n\n${HEADER_ENTRIES}}\n")
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/generated/generated_textures.hpp" "${HEADER_CONTENT}")

add_custom_target(filament_textures DEPENDS ${ALL_KTX_FILES})
add_dependencies(filament_textures mipgen)

# ==============================================================================

file(GLOB_RECURSE SRC_FILES "src/*.cpp")
add_executable(cardillo_frontend ${SRC_FILES})
add_dependencies(cardillo_frontend filament_materials filament_textures)

target_include_directories(cardillo_frontend PRIVATE
    src
    ${CMAKE_CURRENT_SOURCE_DIR}/../src
    ${CMAKE_CURRENT_SOURCE_DIR}/../examples
    ${CMAKE_CURRENT_BINARY_DIR}/generated
)
target_compile_definitions(cardillo_frontend PRIVATE
    FILAMENT_MATERIAL_PATH="${CMAKE_CURRENT_BINARY_DIR}/generated"
    FILAMENT_TEXTURE_PATH="${CMAKE_CURRENT_BINARY_DIR}/generated/textures"
    FILAMENT_IBL_SKYBOX_PATH="${filament_SOURCE_DIR}/docs/web/assets/redball/pillars_2k/pillars_2k_skybox.ktx"
    FILAMENT_IBL_REFLECTION_PATH="${filament_SOURCE_DIR}/docs/web/assets/redball/pillars_2k/pillars_2k_ibl.ktx"
)

target_link_libraries(cardillo_frontend PRIVATE
    filament
    ktxreader
    imageio
    glfw
    filagui
    OpenGL::GL
    cardillo
)

# Only link X11 on Linux
if(NOT WIN32)
    target_link_libraries(cardillo_frontend PRIVATE X11::X11)
endif()

if (MSVC)
    target_compile_options(cardillo_frontend PRIVATE /W4)
    
    target_compile_definitions(cardillo_frontend PRIVATE 
        $<$<CONFIG:Debug>:_ITERATOR_DEBUG_LEVEL=0>
    )
else()
    target_compile_options(cardillo_frontend PRIVATE -Wall -Wextra -Wpedantic)
endif()