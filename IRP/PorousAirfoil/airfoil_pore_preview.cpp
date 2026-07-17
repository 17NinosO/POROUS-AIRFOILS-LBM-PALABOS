#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "config.hpp"
#include "domain.hpp"
#include "lattice.hpp"
#include "sizing.hpp"
#include "naca_geometry.hpp"

using namespace pore_geometry;

namespace {

struct PoreLookup {
    DomainMask mask;
    double pore_domain_width;
    double pore_domain_height;
    double resolution;  // stored explicitly -- DomainMask itself doesn't keep this

    bool isPoreAt(double x_norm, double y_norm) const {
        double x_local = x_norm * pore_domain_width;
        double y_local = y_norm + pore_domain_height / 2.0;
        if (x_local < 0.0 || x_local >= pore_domain_width) return false;
        if (y_local < 0.0 || y_local >= pore_domain_height) return false;
        int col = static_cast<int>(x_local / resolution);
        int row = static_cast<int>(y_local / resolution);
        col = std::clamp(col, 0, mask.nx() - 1);
        row = std::clamp(row, 0, mask.ny() - 1);
        return mask.at(row, col);
    }
};

GeometryConfig baselinePoreConfig() {
    GeometryConfig config;
    config.domain_width = 1.0;
    config.domain_height = 0.16;
    config.buffer_fraction = 0.03;
    config.vertical_margin_fraction = 0.15;

    config.lattice_type = LatticeType::Square;
    config.n_cols = 8;
    config.n_rows = 3;   // multiple rows -> grid of pores through the body
    config.disorder = 0.0;

    config.size_distribution = SizeDistribution::Uniform;
    config.r_mean = 0.012;
    config.r_std = 0;    // effectively zero -- avoids UB in the RNG at
                              // exactly 0, but every pore comes out identical
    config.r_min = 0.012;    // min == max == mean -> exact uniform size
    config.r_max = 0.012;

    config.throat_enabled = true;
    config.throat_reach_domain_edges = true;  // no domain-edge throats for a tapered airfoil
    config.throat_width_mean = 0.004;
    config.throat_width_std = 0;  // same trick -- every throat identical
    config.throat_width_min = 0.002;

    config.n_harmonics = 0;
    config.boundary_roughness = 0.0;

    config.resolution = 0.0005;
    config.seed = 42;
    return config;
}

void writePGM(const std::vector<unsigned char>& solid, int nx, int ny, const std::string& path) {
    std::ofstream out(path, std::ios::binary);
    out << "P5\n" << nx << " " << ny << "\n255\n";
    for (int row = ny - 1; row >= 0; --row)
        for (int col = 0; col < nx; ++col)
            out.put(solid[static_cast<std::size_t>(row) * nx + col] ? 0 : (char)255);
}

}  // namespace

int main() {
    GeometryConfig poreConfig = baselinePoreConfig();
    poreConfig.validate();
    {
    auto centers = generateLatticeCentres(poreConfig);
    auto radii = assignPoreRadii(poreConfig, centers.size());
    std::cout << "raw pore radii: ";
    for (double r : radii) std::cout << r << " ";
    std::cout << "\n";
}
    auto poreMaskRaw = buildDomainMask(poreConfig);
    PoreLookup poreLookup{std::move(poreMaskRaw), poreConfig.domain_width, poreConfig.domain_height, poreConfig.resolution};

    double x0 = -0.05, x1 = 1.05, y0 = -0.12, y1 = 0.12, resolution = 0.0015;
    int nx = static_cast<int>((x1 - x0) / resolution);
    int ny = static_cast<int>((y1 - y0) / resolution);
    std::vector<unsigned char> solid(static_cast<std::size_t>(nx) * ny, 0);

    std::size_t airfoilCount = 0, poreCount = 0;
    for (int row = 0; row < ny; ++row) {
        double y_norm = y0 + row * resolution;
        for (int col = 0; col < nx; ++col) {
            double x_norm = x0 + col * resolution;
            bool inAirfoil = isInsideAirfoil(x_norm, y_norm);
            bool inPore = inAirfoil && poreLookup.isPoreAt(x_norm, y_norm);
            if (inAirfoil) ++airfoilCount;
            if (inPore) ++poreCount;
            solid[static_cast<std::size_t>(row) * nx + col] = (inAirfoil && !inPore) ? 1 : 0;
        }
    }

    writePGM(solid, nx, ny, "airfoil_pores.pgm");
    std::cout << "pore fraction of airfoil interior: " << double(poreCount) / double(airfoilCount) << "\n";
    std::cout << "saved airfoil_pores.pgm (" << nx << "x" << ny << ")\n";
    return 0;
}