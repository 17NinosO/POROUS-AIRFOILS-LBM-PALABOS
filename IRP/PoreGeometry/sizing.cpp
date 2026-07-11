// sizing.cpp
#include "sizing.hpp"
#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

namespace pore_geometry {

std::vector<double> assignPoreRadii(const GeometryConfig& config, std::size_t n_pores) {
    std::mt19937_64 rng(config.seed + 1);
    std::vector<double> radii(n_pores);

    switch (config.size_distribution) {
        case SizeDistribution::Uniform: {
            std::uniform_real_distribution<double> dist(config.r_min, config.r_max);
            for (auto& r : radii) r = dist(rng);
            break;
        }
        case SizeDistribution::Gaussian: {
            std::normal_distribution<double> dist(config.r_mean, config.r_std);
            for (auto& r : radii) r = std::clamp(dist(rng), config.r_min, config.r_max);
            break;
        }
        case SizeDistribution::LogNormal: {
            double variance = config.r_std * config.r_std;
            double mu = std::log(config.r_mean * config.r_mean /
                                  std::sqrt(variance + config.r_mean * config.r_mean));
            double sigma = std::sqrt(std::log(1.0 + variance / (config.r_mean * config.r_mean)));
            std::lognormal_distribution<double> dist(mu, sigma);
            for (auto& r : radii) r = std::clamp(dist(rng), config.r_min, config.r_max);
            break;
        }
        default:
            throw std::invalid_argument("Unknown size distribution");
    }
    return radii;
}

}  // namespace pore_geometry