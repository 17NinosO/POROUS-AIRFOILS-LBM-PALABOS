# PorousChannel/CMakeLists.txt
cmake_minimum_required(VERSION 3.14)
project(porous_channel CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(MPI REQUIRED)
add_definitions(-DPLB_MPI_PARALLEL)

set(PLB_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../..)   # PorousChannel/../.. = palabos/
set(PLB_SRC  ${PLB_ROOT}/src)
set(PLB_EXT  ${PLB_ROOT}/externalLibraries)

file(GLOB_RECURSE PLB_SOURCES   ${PLB_SRC}/*.cpp)
file(GLOB TINYXML_SOURCES ${PLB_EXT}/tinyxml/*.cpp)

include_directories(${PLB_SRC})
include_directories(${PLB_EXT})
include_directories(${MPI_CXX_INCLUDE_DIRS})

add_subdirectory(../PoreGeometry pore_geometry_build)

add_executable(porous_channel
    porous_channel.cpp
    ${PLB_SOURCES}
    ${TINYXML_SOURCES}
)

target_link_libraries(porous_channel PRIVATE ${MPI_CXX_LIBRARIES} m pore_geometry)
target_compile_options(porous_channel PRIVATE -O3 -Wall)