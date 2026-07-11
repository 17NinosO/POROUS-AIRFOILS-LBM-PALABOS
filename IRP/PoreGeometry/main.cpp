// main.cpp
#include <iostream>
#include <algorithm>
#include <numeric>
#include "config.hpp"
#include "lattice.hpp"
#include "sizing.hpp"

int main() {
    pore_geometry::GeometryConfig config;
    config.disorder = 0.0;

    auto centres = pore_geometry::generateLatticeCentres(config);
    auto radii = pore_geometry::assignPoreRadii(config, centres.size());

    double min_r = *std::min_element(radii.begin(), radii.end());
    double max_r = *std::max_element(radii.begin(), radii.end());
    double mean_r = std::accumulate(radii.begin(), radii.end(), 0.0) / radii.size();

    std::cout << "generated " << radii.size() << " radii\n";
    std::cout << "min=" << min_r << " mean=" << mean_r << " max=" << max_r << "\n";
    std::cout << "config range was [" << config.r_min << ", " << config.r_max << "]\n";
    return 0;
}