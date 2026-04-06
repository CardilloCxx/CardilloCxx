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

# If QOCO was installed by top-level ExternalProject targets, use CPU headers
# and runtime-load either CPU/CUDA shared library from their install dirs.
if(NOT DEFINED QOCO_INSTALL_DIR_CPU)
    set(QOCO_INSTALL_DIR_CPU "${CMAKE_BINARY_DIR}/qoco-install-cpu")
endif()
if(NOT DEFINED QOCO_INSTALL_DIR_CUDA)
    set(QOCO_INSTALL_DIR_CUDA "${CMAKE_BINARY_DIR}/qoco-install-cuda")
endif()

target_include_directories(${PROJECT_NAME} PUBLIC "${QOCO_INSTALL_DIR_CPU}/include")
target_include_directories(${PROJECT_NAME} PUBLIC "${QOCO_INSTALL_DIR_CPU}/include/qdldl")
target_include_directories(${PROJECT_NAME} PUBLIC "${CMAKE_BINARY_DIR}/_deps/qoco-cpu/src/qoco_project_cpu/include")

target_compile_definitions(${PROJECT_NAME} PUBLIC HAVE_QOCO=1)
target_compile_definitions(${PROJECT_NAME} PUBLIC QOCO_CPU_LIB_PATH="${QOCO_INSTALL_DIR_CPU}/lib/libqoco.so")

if(EXISTS "${QOCO_INSTALL_DIR_CUDA}/lib/libqoco.so")
    target_compile_definitions(${PROJECT_NAME} PUBLIC QOCO_CUDA_LIB_PATH="${QOCO_INSTALL_DIR_CUDA}/lib/libqoco.so")
endif()

# Optional: enable EXR loading via tinyexr (header-only). Define CARDILLO_USE_TINYEXR to turn it on if header is provided.
# Users can set -DCARDILLO_USE_TINYEXR=ON when configuring to enable.
option(CARDILLO_USE_TINYEXR "Enable EXR heightmap loading via tinyexr" OFF)
if (CARDILLO_USE_TINYEXR)
    target_compile_definitions(${PROJECT_NAME} PUBLIC CARDILLO_USE_TINYEXR=1)
    find_package(ZLIB REQUIRED)
    target_link_libraries(${PROJECT_NAME} PUBLIC ZLIB::ZLIB)
endif()
