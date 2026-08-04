// throats.hpp
#pragma once
#include "lattice.hpp" // Vec2

namespace pore_geometry {

bool inCapsule(Vec2 p, Vec2 p1, Vec2 p2, double half_width);

} //namespace pore_geometry