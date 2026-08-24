// ============================================================================
// generate_conforming_geometry.cpp
//
// Standalone preview of the conforming mini-airfoil lattice. NO Palabos
// dependency -- uses the EXACT SAME parameters as
// NACA0012Palabos_ConformingLattice.cpp's initConformingGeometry(), so
// what you see here is guaranteed to match what the solver actually
// stamps. Check this before spending any HPC time.
//
// Build:
//   g++ -O2 -std=c++14 generate_conforming_geometry.cpp -o generate_conforming_geometry
// Run:
//   ./generate_conforming_geometry
//
// Produces:
//   conforming_centers.csv   -- i,j,cx,cy,theta_deg per mini-airfoil
//   conforming_surfaces.csv  -- upper/lower surface points per mini-airfoil
//                               (same format plot_lattice_geometry.py already reads)
// ============================================================================

#include "AirfoilLatticeGeometry.h"
#include <iostream>

int main() {
    // MUST MATCH initConformingGeometry() in
    // NACA0012Palabos_ConformingLattice.cpp exactly. If you change one,
    // change both, or this preview stops being trustworthy.
    geom::ConformingLatticeParams cp;
    cp.N_chord_big = 2000;          // 10x bigger -- domain in the solver stays UNCHANGED, refine box now hugs actual geometry
    cp.c_mini_frac = 0.0246154;     // solved for ~30 columns
    cp.pitch_s_factor = 1.3;
    cp.pitch_n_factor = 1.0;
    cp.x_margin = 0.02;
    cp.pivot_chords = 0.25;

    geom::LatticeGeometry g = geom::buildConformingLattice(cp);
    geom::plint miniChord = geom::conformingChordUnits(cp);

    auto overlaps = geom::findOverlaps(g, miniChord, cp.pivot_chords);
    if (!overlaps.empty()) {
        std::cout << "WARNING: " << overlaps.size() << " overlapping pair(s).\n";
        for (auto& pr : overlaps) std::cout << "  units[" << pr.first << "] x units[" << pr.second << "]\n";
    } else {
        std::cout << "No overlaps detected.\n";
    }

    std::cout << "Mini-airfoil count: " << g.units.size() << "\n";
    std::cout << "Mini-airfoil chord (coarse lattice units): " << miniChord << "\n";
    std::cout << "Overall/big chord (N_chord_big): " << cp.N_chord_big << "\n";

    geom::writeLatticeCenters(g, "conforming_centers.csv");
    geom::writeLatticeSurfaces(g, miniChord, "conforming_surfaces.csv");
    std::cout << "Wrote conforming_centers.csv and conforming_surfaces.csv\n";

    return 0;
}