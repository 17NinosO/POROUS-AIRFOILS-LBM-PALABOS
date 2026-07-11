// lattice.hpp
#pragma once
#include <array>
#include <vector>
#include "config.hpp"

namespace pore_geometry {

struct Vec2 { double x = 0.0, y = 0.0; };

// Region where the pore centres may be placed: (x0, x1, y0, y1)
std::array<double, 4> porousWindowBounds(const GeometryConfig& config);

// Row-major: index = row * n_cols + col
//this flat list back into a grid to find each pore neighbour
std::vector<Vec2> generateLatticeCentres(const GeometryConfig& config);

} //namespace pore_geometry