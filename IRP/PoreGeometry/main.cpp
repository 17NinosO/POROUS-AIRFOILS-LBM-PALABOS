// main.cpp
#include <iostream>
#include "config.hpp"
#include "domain.hpp"
#include "pgm_writer.hpp"

int main() {
    pore_geometry::GeometryConfig config;
    auto mask = pore_geometry::buildDomainMask(config);
    pore_geometry::writePGM(mask, "baseline_geometry.pgm");
    std::cout << "saved baseline_geometry.pgm (" << mask.nx() << "x" << mask.ny()
              << ", porosity=" << mask.porosity() << ")\n";
    return 0;
}