// config.hpp
#pragma once
#include <stdexcept>

namespace pore_geometry {

    enum class LatticeType { Square, Hex };
    enum class SizeDistribution { Uniform, Gaussian, LogNormal };

struct GeometryConfig {
    // Domain
    double domain_width = 400.0;
    double domain_height = 180.0;
    double buffer_fraction = 0.10;
    double vertical_margin_fraction = 0.10;

    // Lattice pore centre placement
    LatticeType lattice_type = LatticeType::Square;
    int n_cols = 8;
    int n_rows = 5;
    double disorder = 0.0;

    // Pore size distribution
    SizeDistribution size_distribution = SizeDistribution::Uniform;
    double r_mean = 10.0;
    double r_std = 1.0;
    double r_min = 8.0;
    double r_max = 12.0;

    // Throats
    bool throat_enabled = true;
    bool throat_reach_domain_edges = true;  // rectangular-block behaviour by
                                             // default; set false for geometries
                                             // (e.g. an airfoil) where the domain
                                             // boundary isn't a flat wall to open
    double throat_width_mean = 3.0;
    double throat_width_std = 0.5;
    double throat_width_min = 1.5;

    // Organic Boundary
    int n_harmonics = 4;
    double boundary_roughness = 1.0;

    unsigned int seed = 42;
    double resolution = 0.5;

    void validate() const{
        if (disorder < 0.0 || disorder > 1.0)
            throw std::invalid_argument("Disorder must be in [0, 1]");
        if (r_min <= 0.0 || r_max < r_min)
            throw std::invalid_argument("Require 0 < r_min <= r_max");
        // more checks will be added with more parameters
    }
};

}