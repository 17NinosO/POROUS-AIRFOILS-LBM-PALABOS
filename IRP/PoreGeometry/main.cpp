// main.cpp
#include <iostream>
#include "throats.hpp"

int main() {
    pore_geometry::Vec2 p1{0.0, 0.0};
    pore_geometry::Vec2 p2{10.0, 0.0};
    double half_width = 2.0;

    struct Check { pore_geometry::Vec2 p; bool expected; const char* label; };
    Check checks[] = {
        {{5.0, 0.0},  true,  "midpoint of segment"},
        {{5.0, 1.9},  true,  "just inside, above midpoint"},
        {{5.0, 2.1},  false, "just outside, above midpoint"},
        {{0.0, 0.0},  true,  "exactly at p1"},
        {{-0.5, 0.0}, true,  "past p1, but inside its rounded cap"},
        {{-2.0, 0.0}, true,  "exactly on the cap boundary (distance == half_width)"},
        {{-2.1, 0.0}, false, "just outside the cap"},
        {{-3.0, 0.0}, false, "well past p1, outside the capsule entirely"},
        {{13.0, 0.0}, false, "well past p2, outside the capsule entirely"},
    };

    for (auto& c : checks) {
        bool got = pore_geometry::inCapsule(c.p, p1, p2, half_width);
        std::cout << c.label << ": got=" << got << " expected=" << c.expected
                   << (got == c.expected ? "  OK" : "  MISMATCH") << "\n";
    }
    return 0;
}