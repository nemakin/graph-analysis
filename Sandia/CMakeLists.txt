cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project(sandia_spla CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add SPLA subdirectory
add_subdirectory(../spla spla_build)

# Create executable
add_executable(sandia_spla sandia_spla.cpp)

# Include directories for the target
target_include_directories(sandia_spla PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../spla/include)

# Link against SPLA
target_link_libraries(sandia_spla PRIVATE spla)

# Set optimization flags
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(sandia_spla PRIVATE -O3)
endif()

# Set output directory
set_target_properties(sandia_spla PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
