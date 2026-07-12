# PoreGeometry/CMakeLists.txt
cmake_minimum_required(VERSION 3.14)
project(pore_geometry CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(pore_geometry STATIC
    lattice.cpp
    sizing.cpp
    throats.cpp
    domain.cpp
)
target_include_directories(pore_geometry PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})