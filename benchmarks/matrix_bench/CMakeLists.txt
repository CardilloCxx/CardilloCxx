cmake_minimum_required(VERSION 3.16)

add_executable(matrix_bench
    src/main.cpp
)

# Headers only in this subproject

target_include_directories(matrix_bench PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
)

# Link Eigen for fair comparison and reuse types

target_link_libraries(matrix_bench PRIVATE
    Eigen3::Eigen
)

# Try to find BLAS (CBLAS)
find_package(BLAS)
if(BLAS_FOUND)
    target_compile_definitions(matrix_bench PRIVATE BENCH_HAVE_BLAS=1)
    target_link_libraries(matrix_bench PRIVATE ${BLAS_LIBRARIES})
endif()

# Optionally find CUDA toolkit for cuBLAS; we won't add CUDA sources here, just link
# Try to locate CUDA toolkit for potential cuBLAS backend (host-only, no nvcc needed)
find_package(CUDAToolkit QUIET)
if (CUDAToolkit_FOUND)
    message(STATUS "CUDAToolkit found: enabling cuBLAS backend")
    target_link_libraries(matrix_bench PRIVATE CUDA::cublas CUDA::cudart)
    target_compile_definitions(matrix_bench PRIVATE BENCH_HAVE_CUBLAS)
    # cuSPARSE for sparse GPU SpMV
    if (TARGET CUDA::cusparse)
        message(STATUS "CUDA cuSPARSE found: enabling cuSPARSE backend")
        target_link_libraries(matrix_bench PRIVATE CUDA::cusparse)
        target_compile_definitions(matrix_bench PRIVATE BENCH_HAVE_CUSPARSE)
    endif()
endif()

set_target_properties(matrix_bench PROPERTIES OUTPUT_NAME matrix_bench)

# Optional: gprof profiling support for this target
option(BENCH_PROFILING "Enable -pg profiling for matrix_bench" OFF)
if (BENCH_PROFILING)
    message(STATUS "Enabling gprof profiling (-pg) for matrix_bench")
    target_compile_options(matrix_bench PRIVATE -pg)
    target_link_options(matrix_bench PRIVATE -pg)
endif()
