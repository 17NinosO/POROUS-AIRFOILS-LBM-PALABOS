// ============================================================================
// generate_lattice_geometry.cpp
//
// Standalone geometry generator. NO Palabos dependency -- builds on the
// pure math in AirfoilLatticeGeometry.h and writes out the lattice so you
// can check/plot it (e.g. with plot_validation.py-style scripts) before
// it ever touches Crescent2 or the solver.
//
// Build:
//   g++ -O2 -std=c++14 generate_lattice_geometry.cpp -o generate_lattice_geometry
// Run:
//   ./generate_lattice_geometry
//
// Produces:
//   lattice_centers.csv   -- one row per airfoil: i,j,cx,cy,theta_deg
//   lattice_surfaces.csv  -- upper/lower surface points per airfoil,
//                            same convention as writeAirfoilGeometry()
// ============================================================================

#include "AirfoilLatticeGeometry.h"
#include <iostream>
#include <cmath>

int main() {
    geom::LatticeParams p;

    // --- grid shape ---
    p.nx = 6;
    p.ny = 4;
    p.N_chord = 200;

    // --- spacing ---
    // X: 1/4-chord GAP between one airfoil's TE and the next's LE.
    //    pitch = chord (1.0) + gap (0.25)
    p.pitch_x_chords = 1.25;
    // Y: 1/2-chord GAP between rows, measured against NACA0012's max
    //    thickness (0.12c) since that's the airfoil's own extent in y.
    //    pitch = max_thickness (0.12) + gap (0.5)
    p.pitch_y_chords = 0.62;

    // --- rotation pivot: quarter-chord, standard aerodynamic convention ---
    p.pivot_chords = 0.25;

    // --- boundary clearance, scaled to the WHOLE ARRAY's own footprint
    // (treating the lattice as one bluff body) rather than a single
    // chord. These multipliers (8 / 27 / 8) are the EXACT proportions
    // from the original single-airfoil solver (x_foil=8*N_chord,
    // Lx=36*N_chord, Ly=16*N_chord -> inlet=8c, outlet=27c, wall=8c each
    // side), now applied to the array's own width/height instead of one
    // chord. This intentionally produces a MUCH larger absolute domain
    // than the earlier 3/3/1.5 version -- HPC-scale only, not for local
    // smoke testing at full nx x ny. ---
    // --- boundary clearance, scaled to the WHOLE ARRAY's own footprint
    // (treating the lattice as one bluff body). Keeps the original
    // solver's asymmetry (inlet closer, outlet further for wake) at
    // multipliers that stay compute-feasible for the array's own size. ---
    geom::applyBodyScaledMargins(p, /*inlet_widths=*/2.0, /*outlet_widths=*/4.0, /*wall_heights=*/2.0);

    // --- angle field: replace with a lookup table / optimiser output
    // whenever you're ready to move past this placeholder. This example
    // is a smooth spatial variation so the "dynamic porous media"
    // character is visible immediately. ---
    p.theta_deg_fn = [](int i, int j) {
        return 15.0 * std::sin(0.6 * j) * std::cos(0.4 * i);
    };

    geom::LatticeGeometry g = geom::buildLattice(p);

    std::vector<std::pair<int,int>> overlaps = geom::findOverlaps(g, p.N_chord, p.pivot_chords);
    if (!overlaps.empty()) {
        std::cout << "WARNING: " << overlaps.size()
                  << " airfoil pair(s) geometrically overlap at this pitch/angle "
                  << "combination.\n";
        for (auto& pr : overlaps) {
            std::cout << "  units[" << pr.first << "] x units[" << pr.second << "]\n";
        }
    } else {
        std::cout << "No overlaps detected.\n";
    }

    std::cout << "Resolved margins (chords): inlet=" << p.margin_inlet_chords
              << " outlet=" << p.margin_outlet_chords
              << " wall=" << p.margin_wall_chords << "\n";
    std::cout << "Lattice: " << p.nx << " x " << p.ny << " airfoils\n";
    std::cout << "Chord N_chord = " << p.N_chord << " lattice units\n";
    std::cout << "Derived domain (coarse lattice units): Lx = " << g.Lx
              << ", Ly = " << g.Ly << "\n";
    std::cout << "Derived domain (chords): Lx = "
              << static_cast<double>(g.Lx) / p.N_chord << "C, Ly = "
              << static_cast<double>(g.Ly) / p.N_chord << "C\n";

    geom::writeLatticeCenters(g, "lattice_centers.csv");
    geom::writeLatticeSurfaces(g, p.N_chord, "lattice_surfaces.csv");

    std::cout << "Wrote lattice_centers.csv and lattice_surfaces.csv\n";
    return 0;
}