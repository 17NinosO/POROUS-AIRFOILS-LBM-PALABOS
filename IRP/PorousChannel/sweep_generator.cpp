// sweep_generator.cpp
//
// Generates 10 symmetric geometries per parameter axis (throat width,
// pore size, pore packing density, pore shape, disorder/"uniformity"),
// holding every other parameter at a fixed symmetric baseline. Only the
// uniformity sweep intentionally breaks symmetry.
//
// Outputs:
//   sweep_output/<axis>_<index>.pgm   -- one image per geometry
//   sweep_output/sweep_manifest.csv   -- axis, index, swept value, porosity

#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include "config.hpp"
#include "domain.hpp"
#include "pgm_writer.hpp"

using namespace pore_geometry;
namespace fs = std::filesystem;

namespace {

GeometryConfig baselineConfig() {
    GeometryConfig config;
    config.domain_width = 400.0;
    config.domain_height = 180.0;
    config.lattice_type = LatticeType::Square;
    config.n_cols = 8;
    config.n_rows = 5;
    config.disorder = 0.0;              // symmetric baseline

    config.size_distribution = SizeDistribution::Uniform;
    config.r_mean = 10.0;
    config.r_std = 0.1;                 // near-zero spread -> symmetric sizes
    config.r_min = 8.0;
    config.r_max = 12.0;

    config.throat_width_mean = 3.0;
    config.throat_width_std = 0.1;      // near-zero spread -> symmetric throats
    config.throat_width_min = 0.5;

    config.n_harmonics = 4;
    config.boundary_roughness = 0.0;    // baseline = smooth circular pores

    config.resolution = 0.5;
    config.seed = 42;
    return config;
}

std::vector<double> linspace(double lo, double hi, int n) {
    std::vector<double> v(n);
    for (int i = 0; i < n; ++i)
        v[i] = lo + (hi - lo) * static_cast<double>(i) / static_cast<double>(n - 1);
    return v;
}

std::string indexedName(const std::string& axis, int i) {
    std::string idx = (i < 10 ? "0" : "") + std::to_string(i);
    return "sweep_output/" + axis + "_" + idx + ".pgm";
}

void runSweep(const std::string& axisName, const std::vector<double>& values,
              std::ofstream& manifest,
              const std::function<void(GeometryConfig&, double)>& apply) {
    for (std::size_t i = 0; i < values.size(); ++i) {
        GeometryConfig config = baselineConfig();
        apply(config, values[i]);
        config.validate();

        auto mask = buildDomainMask(config);
        std::string filename = indexedName(axisName, static_cast<int>(i));
        writePGM(mask, filename);

        manifest << axisName << "," << i << "," << values[i] << "," << mask.porosity() << "\n";
        std::cout << axisName << " [" << i << "] value=" << values[i]
                   << " porosity=" << mask.porosity() << " -> " << filename << "\n";
    }
}

}  // namespace

int main() {
    fs::create_directories("sweep_output");
    std::ofstream manifest("sweep_output/sweep_manifest.csv");
    manifest << "axis,index,value,porosity\n";

    runSweep("throat_width", linspace(1.0, 8.0, 10), manifest,
             [](GeometryConfig& c, double v) { c.throat_width_mean = v; });

    runSweep("pore_size", linspace(6.0, 16.0, 10), manifest,
             [](GeometryConfig& c, double v) {
                 c.r_mean = v;
                 c.r_min = v - 2.0;
                 c.r_max = v + 2.0;
             });

    // Pore "distribution" interpreted as packing density (n_cols is an
    // integer, so this sweeps 4..13 pores per row, not a continuous value).
    for (int i = 0; i < 10; ++i) {
        GeometryConfig config = baselineConfig();
        config.n_cols = 4 + i;
        config.validate();
        auto mask = buildDomainMask(config);
        std::string filename = indexedName("pore_distribution", i);
        writePGM(mask, filename);
        manifest << "pore_distribution," << i << "," << config.n_cols << "," << mask.porosity() << "\n";
        std::cout << "pore_distribution [" << i << "] n_cols=" << config.n_cols
                  << " porosity=" << mask.porosity() << " -> " << filename << "\n";
    }

    runSweep("pore_shape", linspace(0.0, 6.0, 10), manifest,
             [](GeometryConfig& c, double v) { c.boundary_roughness = v; });

    // The one axis that intentionally breaks symmetry, by design.
    runSweep("uniformity", linspace(0.0, 1.0, 10), manifest,
             [](GeometryConfig& c, double v) { c.disorder = v; });

    std::cout << "\nDone. 50 geometries written to sweep_output/, manifest at sweep_output/sweep_manifest.csv\n";
    return 0;
}