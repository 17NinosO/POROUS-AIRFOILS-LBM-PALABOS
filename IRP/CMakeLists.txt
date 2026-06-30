cmake_minimum_required(VERSION 3.0)
project(NACA0012 CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# MPI
find_package(MPI REQUIRED)
add_definitions(-DPLB_MPI_PARALLEL)

# Palabos paths
set(PLB_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/..)
set(PLB_SRC  ${PLB_ROOT}/src)
set(PLB_EXT  ${PLB_ROOT}/externalLibraries)   # ← new

# Collect sources
file(GLOB_RECURSE PLB_SOURCES   ${PLB_SRC}/*.cpp)
file(GLOB TINYXML_SOURCES ${PLB_EXT}/tinyxml/*.cpp)

# Include directories
include_directories(${PLB_SRC})
include_directories(${PLB_EXT})                     # ← new
include_directories(${MPI_CXX_INCLUDE_DIRS})

# Executable
add_executable(NACA0012
    NACA0012Palabos.cpp
    ${PLB_SOURCES}
    ${TINYXML_SOURCES}
)

target_link_libraries(NACA0012 PRIVATE ${MPI_CXX_LIBRARIES} m)
target_compile_options(NACA0012 PRIVATE -O3 -Wall)