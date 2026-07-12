// domain.cpp
#include "domain.hpp"
#include <algorithm>
#include <cmath>
#include <random>
#include "lattice.hpp"
#include "sizing.hpp"
#include "pore_shape.hpp"
#include "throats.hpp"

namespace pore_geometry {

namespace {

double throatHalfWidth(const GeometryConfig& config, double r1, double r2, std::mt19937_64& rng) {
    std::normal_distribution<double> dist(config.throat_width_mean, config.throat_width_std);
    double width = std::clamp(dist(rng), config.throat_width_min, 2.0 * std::min(r1, r2));
    return width / 2.0;
}

}  // namespace

DomainMask buildDomainMask(const GeometryConfig& config) {
    config.validate();

    int nx = std::max(int(std::round(config.domain_width / config.resolution)), 2);
    int ny = std::max(int(std::round(config.domain_height / config.resolution)), 2);
    DomainMask mask(nx, ny);

    auto centers = generateLatticeCentres(config);          // <- match YOUR spelling
    auto radii = assignPoreRadii(config, centers.size());

    std::mt19937_64 shape_rng(config.seed + 2);
    std::vector<PoreShape> pores;
    pores.reserve(centers.size());
    for (std::size_t k = 0; k < centers.size(); ++k)
        pores.push_back(makeOrganicPore(centers[k], radii[k], config.n_harmonics,
                                         config.boundary_roughness, shape_rng));

    // carve out pores
    for (int row = 0; row < ny; ++row) {
        double y = row * config.resolution;
        for (int col = 0; col < nx; ++col) {
            double x = col * config.resolution;
            for (auto& pore : pores) {
                if (pore.contains({x, y})) { mask.set(row, col, true); break; }
            }
        }
    }

    // carve out throats between grid neighbours only
    if (config.throat_enabled) {
        std::mt19937_64 throat_rng(config.seed + 3);
        auto index = [&](int row, int col) { return std::size_t(row) * config.n_cols + col; };

        for (int j = 0; j < config.n_rows; ++j) {
            for (int i = 0; i < config.n_cols; ++i) {
                if (i + 1 < config.n_cols) {
                    auto a = index(j, i), b = index(j, i + 1);
                    double half_w = throatHalfWidth(config, radii[a], radii[b], throat_rng);
                    for (int row = 0; row < ny; ++row)
                        for (int col = 0; col < nx; ++col)
                            if (inCapsule({col * config.resolution, row * config.resolution},
                                          centers[a], centers[b], half_w))
                                mask.set(row, col, true);
                }
                if (j + 1 < config.n_rows) {
                    auto a = index(j, i), b = index(j + 1, i);
                    double half_w = throatHalfWidth(config, radii[a], radii[b], throat_rng);
                    for (int row = 0; row < ny; ++row)
                        for (int col = 0; col < nx; ++col)
                            if (inCapsule({col * config.resolution, row * config.resolution},
                                          centers[a], centers[b], half_w))
                                mask.set(row, col, true);
                }
            }
        }
    }

    return mask;
}

}  // namespace pore_geometry