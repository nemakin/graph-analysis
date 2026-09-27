cmake_minimum_required(VERSION 3.23)
project(tc LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SPLA_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(SPLA_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(SPLA_BUILD_OPENCL ON CACHE BOOL "" FORCE)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../spla ${CMAKE_CURRENT_BINARY_DIR}/spla)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../GraphBLAS ${CMAKE_CURRENT_BINARY_DIR}/GraphBLAS)

add_library(tc_algorithms
    src/graphblas/tc.cpp
    src/graphblas/utils.cpp
    src/spla/tc.cpp
    src/spla/utils.cpp
)

target_link_libraries(tc_algorithms PUBLIC
        spla
        GraphBLAS
)

target_include_directories(tc_algorithms PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_executable(tc_bench tc_bench.cpp)
target_link_libraries(tc_bench PRIVATE tc_algorithms)

add_executable(tc_profile tc_profile.cpp)
target_link_libraries(tc_profile PRIVATE tc_algorithms)

include(CTest)
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
