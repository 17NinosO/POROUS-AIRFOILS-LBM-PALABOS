cmake_minimum_required(VERSION 3.14)
project(porous_channel CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(../PoreGeometry pore_geometry_build)

add_executable(bridge_test bridge_test.cpp)
target_link_libraries(bridge_test PRIVATE pore_geometry)