// main.cpp
#include <iostream>
#include "config.hpp"
#include "domain.hpp"

int main() {
    pore_geometry::GeometryConfig config;
    auto mask = pore_geometry::buildDomainMask(config);
    std::cout << "grid: " << mask.nx() << " x " << mask.ny() << "\n";
    std::cout << "porosity: " << mask.porosity() << "\n";
    return 0;
}