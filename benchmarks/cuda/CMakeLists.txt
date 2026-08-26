cmake_minimum_required(VERSION 3.20)
project(xfemm_cuda_benchmark LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(xfemm-system-inspect
    inspect_system.cpp
    xfemm_system.cpp
)
target_compile_options(xfemm-system-inspect PRIVATE -Wall -Wextra -Wpedantic)

add_executable(xfemm-replace-ans-solution replace_ans_solution.cpp)
target_compile_options(xfemm-replace-ans-solution PRIVATE -Wall -Wextra -Wpedantic)

add_executable(xfemm-age-topology-enumerate
    age_topology_enumerate.cpp
    xfemm_system.cpp
)
target_compile_options(xfemm-age-topology-enumerate PRIVATE -Wall -Wextra -Wpedantic)

option(XFEMM_CUDA_BENCHMARK "Build the CUDA solver benchmark" ON)
if(XFEMM_CUDA_BENCHMARK)
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES native)
    endif()
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    add_executable(xfemm-cuda-benchmark
        benchmark.cu
        xfemm_system.cpp
    )
    set_target_properties(xfemm-cuda-benchmark PROPERTIES
        CUDA_STANDARD 17
        CUDA_STANDARD_REQUIRED ON
    )
    target_compile_options(xfemm-cuda-benchmark PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic>
    )
    target_link_libraries(xfemm-cuda-benchmark PRIVATE
        CUDA::cudart CUDA::cublas CUDA::cusparse
    )

    add_executable(xfemm-cuda-ssor-benchmark
        ssor_benchmark.cu
        xfemm_system.cpp
    )
    set_target_properties(xfemm-cuda-ssor-benchmark PROPERTIES
        CUDA_STANDARD 17
        CUDA_STANDARD_REQUIRED ON
    )
    target_compile_options(xfemm-cuda-ssor-benchmark PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic>
    )
    target_link_libraries(xfemm-cuda-ssor-benchmark PRIVATE
        CUDA::cudart CUDA::cublas CUDA::cusparse
    )

    set(CUDSS_ROOT "" CACHE PATH "Optional cuDSS installation prefix")
    find_path(CUDSS_INCLUDE_DIR cudss.h HINTS "${CUDSS_ROOT}"
        PATHS /usr/include/libcudss/13
        PATH_SUFFIXES include include/libcudss/13)
    find_library(CUDSS_LIBRARY NAMES cudss libcudss.so.0
        HINTS "${CUDSS_ROOT}"
        PATHS /usr/lib/x86_64-linux-gnu/libcudss/13
        PATH_SUFFIXES lib lib64 lib/x86_64-linux-gnu/libcudss/13)
    if(CUDSS_INCLUDE_DIR AND CUDSS_LIBRARY)
        target_sources(xfemm-cuda-benchmark PRIVATE cudss_solver.cpp)
        target_include_directories(xfemm-cuda-benchmark PRIVATE "${CUDSS_INCLUDE_DIR}")
        target_link_libraries(xfemm-cuda-benchmark PRIVATE "${CUDSS_LIBRARY}")
        target_compile_definitions(xfemm-cuda-benchmark PRIVATE XFEMM_HAVE_CUDSS=1)

        add_executable(xfemm-cuda-age-topology-benchmark
            age_topology_benchmark.cpp
            cudss_solver.cpp
            xfemm_system.cpp
        )
        target_include_directories(xfemm-cuda-age-topology-benchmark PRIVATE
            "${CUDSS_INCLUDE_DIR}"
            "${CUDAToolkit_INCLUDE_DIRS}"
        )
        target_compile_options(xfemm-cuda-age-topology-benchmark PRIVATE
            -Wall -Wextra -Wpedantic)
        target_link_libraries(xfemm-cuda-age-topology-benchmark PRIVATE
            CUDA::cudart "${CUDSS_LIBRARY}")
        message(STATUS "cuDSS benchmark enabled: ${CUDSS_LIBRARY}")
    else()
        message(STATUS "cuDSS not found; building cuSPARSE PCG only")
    endif()
endif()
