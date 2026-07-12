// throats.cpp
#include "throats.hpp"
#include <cmath>
#include <algorithm>

namespace pore_geometry {

bool inCapsule(Vec2 p, Vec2 p1, Vec2 p2, double half_width) {
    Vec2 seg{p2.x - p1.x, p2.y - p1.y};
    double seg_len_sq = seg.x * seg.x + seg.y * seg.y;

    if (seg_len_sq < 1e-12) {

        return std::hypot(p.x - p1.x, p.y - p1.y) <= half_width;
    }

double t = ((p.x - p1.x) * seg.x + (p.y - p1.y) * seg.y) / seg_len_sq;
    t = std::clamp(t, 0.0, 1.0);

    double closest_x = p1.x + t * seg.x;
    double closest_y = p1.y + t * seg.y;
    return std::hypot(p.x - closest_x, p.y - closest_y) <= half_width;
}

} //namespace pore_geometry