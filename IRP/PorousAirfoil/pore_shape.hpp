#pragma once
#include <cmath>
#include <random>
#include <vector>
#include "lattice.hpp" //for Vec2

namespace pore_geometry {

class PoreShape {
public:
    PoreShape(Vec2 centre, double base_radius,
              std::vector<double> amplitudes, std::vector<double> phases)
            : centre_(centre), base_radius_(base_radius),
              amplitudes_(std::move(amplitudes)), phases_(std::move(phases)) {}

    double radiusAt(double theta) const {
        double r = base_radius_;
        for (std::size_t k = 0; k < amplitudes_.size(); ++k) {
            r += amplitudes_[k] * std::sin(static_cast<double>(k + 1) * theta + phases_[k]);
        }
        return r;
    }

    bool contains(Vec2 p) const {
        Vec2 d{p.x - centre_.x, p.y - centre_.y};
        double theta = std::atan2(d.y, d.x);
        double dist = std::hypot(d.x, d.y);
        return dist <= radiusAt(theta);
    }

private:
    Vec2 centre_;
    double base_radius_;
    std::vector<double> amplitudes_;
    std::vector<double> phases_;
};

// Builds a PoreShape with random amplitudes and phases for the harmonics.
// Takes a random number generator, a base radius, and the number of harmonics to use.

template <class Rng>
PoreShape makeOrganicPore(Vec2 centre, double base_radius, int n_harmonics,
                          double boundary_roughness, Rng& rng) {
    std::vector<double> amplitudes(n_harmonics), phases(n_harmonics);
    std::uniform_real_distribution<double> amp_dist(-1.0, 1.0);
    std::uniform_real_distribution<double> phase_dist(0.0, 2.0 * M_PI);

    for (int k = 0; k < n_harmonics; ++k) {
        double max_amp = boundary_roughness / (k + 1);
        amplitudes[k] = amp_dist(rng) * max_amp;
        phases[k] = phase_dist(rng);
    }
    return PoreShape(centre, base_radius, std::move(amplitudes), std::move(phases));
}

} //namespace pore_geometry