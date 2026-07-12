// PorousChannel/bridge_test.cpp
#include <cmath>
#include <iostream>
#include "config.hpp"
#include "domain.hpp"
#include "pore_mask_domain.hpp"

int main() {
    pore_geometry::GeometryConfig config;
    config.disorder = 0.0;
    config.r_mean = 10.0;
    config.r_std = 1.0;
    config.r_min = 8.0;
    config.r_max = 12.0;
    config.throat_width_mean = 3.0;
    config.throat_width_std = 0.5;
    config.throat_width_min = 1.0;
    config.boundary_roughness = 1.0;

    auto mask = pore_geometry::buildDomainMask(config);
    PoreMaskDomain domain(mask);

    std::size_t solid_count = 0;
    for (int iY = 0; iY < domain.ny(); ++iY)
        for (int iX = 0; iX < domain.nx(); ++iX)
            if (domain.isSolid(iX, iY)) ++solid_count;

    double wrapper_porosity = 1.0 - double(solid_count) / (double(domain.nx()) * domain.ny());

    std::cout << "mask.porosity()  = " << mask.porosity() << "\n";
    std::cout << "wrapper porosity = " << wrapper_porosity << "\n";
    std::cout << (std::abs(mask.porosity() - wrapper_porosity) < 1e-9 ? "MATCH\n" : "MISMATCH -- indexing bug, stop here\n");

    return 0;
}