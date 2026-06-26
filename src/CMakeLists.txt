file(GLOB_RECURSE CARDILLO_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
)

add_library(${PROJECT_NAME} STATIC ${CARDILLO_SOURCES})

target_include_directories(${PROJECT_NAME} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${entt_SOURCE_DIR}/src
)


target_link_libraries(${PROJECT_NAME}
    PUBLIC
        Eigen3::Eigen
        dl
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
target_compile_definitions(${PROJECT_NAME} PUBLIC QOCO_CUDA_LIB_PATH="${QOCO_INSTALL_DIR_CUDA}/lib/libqoco.so")

target_compile_definitions(${PROJECT_NAME} PUBLIC HAVE_CLARABEL=1)
if(TARGET libclarabel_c_shared)
    target_link_libraries(${PROJECT_NAME} PUBLIC libclarabel_c_shared)
endif()


