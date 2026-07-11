#include <iostream>
#include "config.hpp"
#include "lattice.hpp"

int main() {
    pore_geometry::GeometryConfig config;
    config.disorder = 0.0;

    auto centres = pore_geometry::generateLatticeCentres(config);

    std::cout << "generated " << centres.size() << " pore centres:\n";
    std::cout << "first row:\n";
    for (int i = 0; i < config.n_cols; ++i) {
        std::cout << " (" << centres[i].x << ", " << centres[i].y << ")\n";
    }

    return 0;
}