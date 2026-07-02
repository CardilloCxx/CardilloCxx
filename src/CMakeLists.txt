cmake_minimum_required(VERSION 3.16)
project(cardillo LANGUAGES CXX)
set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "Install path prefix" FORCE)
include(FetchContent)

# --- 1. Dependencies ---
set(BUILD_PYTHON_INTERFACE OFF CACHE BOOL "Disable Coal Python bindings" FORCE)
set(COAL_BUILD_PYTHON_INTERFACE OFF CACHE BOOL "Disable Coal Python bindings" FORCE)
set(COAL_BUILD_TESTING OFF CACHE BOOL "Disable Coal testing" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)

FetchContent_Declare(eigen GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git GIT_TAG 3.4.0)
FetchContent_Declare(entt GIT_REPOSITORY https://github.com/skypjack/entt.git GIT_TAG v3.13.2)
FetchContent_Declare(clarabel GIT_REPOSITORY https://github.com/oxfordcontrol/Clarabel.cpp.git GIT_TAG v0.11.1)
FetchContent_Declare(coal GIT_REPOSITORY https://github.com/coal-library/coal.git GIT_TAG v3.0.0)
FetchContent_Declare(qoco GIT_REPOSITORY https://github.com/qoco-org/qoco.git GIT_TAG v0.3.2)
FetchContent_MakeAvailable(eigen coal entt clarabel qoco)

# --- 2. Engine Targets ---
file(GLOB_RECURSE CARDILLO_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
add_library(cardillo SHARED ${CARDILLO_SOURCES})

# --- 3. Link Dependencies ---
target_link_libraries(cardillo PUBLIC 
    Eigen3::Eigen 
    EnTT::EnTT 
    coal::coal
    qoco
    libclarabel_c_shared
)

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