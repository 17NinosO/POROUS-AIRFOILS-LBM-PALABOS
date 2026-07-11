// lattice.cpp
#include "lattice.hpp"
#include <random>

namespace pore_geometry {

std::array<double, 4> porousWindowBounds(const GeometryConfig& config) {
    double x0 = config.domain_width * config.buffer_fraction;
    double x1 = config.domain_width * (1.0 - config.buffer_fraction);
    double y0 = config.domain_height * config.vertical_margin_fraction;
    double y1 = config.domain_height * (1.0 - config.vertical_margin_fraction);
    return {x0, x1, y0, y1};

}

std::vector<Vec2> generateLatticeCentres(const GeometryConfig& config) {
    std::mt19937_64 rng(config.seed);
    std::uniform_real_distribution<double> u(0.0, 1.0);

    auto [x0, x1, y0, y1] = porousWindowBounds(config);
    double spacing_x = (x1 - x0) / (config.n_cols);
    double spacing_y = (y1 - y0) / (config.n_rows);

    std::vector<Vec2> centres;
    centres.reserve(config.n_rows * config.n_cols);

    for (int j = 0; j < config.n_rows; ++j) {
        for (int i = 0; i < config.n_cols; ++i) {
            double x = x0 + (i + 0.5) * spacing_x;
            double y = y0 + (j + 0.5) * spacing_y;
            
            if (config.lattice_type == LatticeType::Hex && (j % 2 == 1))
                x += spacing_x * 0.5; // Offset for hexagonal lattice

            if (config.disorder > 0.0) {
                x += config.disorder * (u(rng) - 0.5) * spacing_x;
                y += config.disorder * (u(rng) - 0.5) * spacing_y;
            }

            centres.push_back({x, y});
        }
    }

    return centres;
}

} //namespace pore_geometry