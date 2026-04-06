file(GLOB_RECURSE CARDILLO_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
)

add_library(${PROJECT_NAME} STATIC ${CARDILLO_SOURCES})

target_include_directories(${PROJECT_NAME} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${entt_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/lib/third_party/tinyexr
)

target_link_libraries(${PROJECT_NAME}
    PUBLIC
        Eigen3::Eigen
        PkgConfig::PETSC
        dl
    MPI::MPI_CXX
)

add_dependencies(${PROJECT_NAME} coal_project)

# Include and link directories after install
target_include_directories(${PROJECT_NAME} PUBLIC
    "${COAL_INSTALL_DIR}/cmeel.prefix/include/"
)

target_link_directories(${PROJECT_NAME} PUBLIC
    "${COAL_INSTALL_DIR}/cmeel.prefix/lib/"
)

target_link_libraries(${PROJECT_NAME} INTERFACE coal)

target_compile_definitions(${PROJECT_NAME} PUBLIC CARDILLO_HAVE_COAL=1)

# If QOCO was installed by the top-level ExternalProject, link against it
if(NOT DEFINED QOCO_INSTALL_DIR)
    set(QOCO_INSTALL_DIR "${CMAKE_BINARY_DIR}/qoco-install")
endif()
# Add QOCO include dir (may be populated by ExternalProject at build time)
target_include_directories(${PROJECT_NAME} PUBLIC "${QOCO_INSTALL_DIR}/include")
target_include_directories(${PROJECT_NAME} PUBLIC "${QOCO_INSTALL_DIR}/include/qdldl")
target_include_directories(${PROJECT_NAME} PUBLIC "${CMAKE_BINARY_DIR}/_deps/qoco/src/qoco_project/include")

message(STATUS "Linking QOCO from ${QOCO_INSTALL_DIR}")
target_link_libraries(${PROJECT_NAME} PUBLIC "${QOCO_INSTALL_DIR}/lib/libqoco.so")
target_compile_definitions(${PROJECT_NAME} PUBLIC HAVE_QOCO=1)

if(DEFINED QOCO_CUDSS_LIBRARY_DIR AND NOT "${QOCO_CUDSS_LIBRARY_DIR}" STREQUAL "")
    target_compile_definitions(${PROJECT_NAME} PUBLIC QOCO_CUDSS_LIBRARY_DIR="${QOCO_CUDSS_LIBRARY_DIR}")
endif()
if(DEFINED QOCO_CUDA_LIBRARY_DIR_RESOLVED AND NOT "${QOCO_CUDA_LIBRARY_DIR_RESOLVED}" STREQUAL "")
    target_compile_definitions(${PROJECT_NAME} PUBLIC QOCO_CUDA_LIBRARY_DIR="${QOCO_CUDA_LIBRARY_DIR_RESOLVED}")
elseif(DEFINED CUDAToolkit_LIBRARY_DIR AND NOT "${CUDAToolkit_LIBRARY_DIR}" STREQUAL "")
    target_compile_definitions(${PROJECT_NAME} PUBLIC QOCO_CUDA_LIBRARY_DIR="${CUDAToolkit_LIBRARY_DIR}")
endif()

# Optional: enable EXR loading via tinyexr (header-only). Define CARDILLO_USE_TINYEXR to turn it on if header is provided.
# Users can set -DCARDILLO_USE_TINYEXR=ON when configuring to enable.
option(CARDILLO_USE_TINYEXR "Enable EXR heightmap loading via tinyexr" OFF)
if (CARDILLO_USE_TINYEXR)
    target_compile_definitions(${PROJECT_NAME} PUBLIC CARDILLO_USE_TINYEXR=1)
    find_package(ZLIB REQUIRED)
    target_link_libraries(${PROJECT_NAME} PUBLIC ZLIB::ZLIB)
endif()
