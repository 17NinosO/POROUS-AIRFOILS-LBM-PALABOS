// ============================================================================
// AirfoilLatticeGeometry.h
//
// Standalone geometry module for a parametrised NxM lattice of NACA0012
// airfoils, each independently rotated by its own angle theta_ij. This is
// the "dynamic porous media" geometry discussed alongside your solid /
// porous single-airfoil work in NACA0012Palabos_Multigrid.cpp.
//
// DELIBERATELY has NO Palabos dependency (no palabos2D.h, no plb::, no
// MultiBlockLattice2D). It only knows about points, chords, and angles.
// This is the separation you asked for: this file compiles and runs on
// its own with a plain C++ compiler, so you can inspect/plot a lattice
// layout in seconds without touching the solver or the HPC queue.
//
// The solver later #includes this header and calls isInsideLattice() from
// inside its own BoxProcessingFunctional2D_S<int> stamping class -- the
// same role isInsideAirfoilReal() plays in NACA0012Palabos_Multigrid.cpp.
// See the comment block at the bottom of this file for that hook.
//
// DESIGN NOTE -- boundary margins are enforced by construction, not by
// manual tuning:
//   You supply margin_inlet_chords, margin_outlet_chords, and
//   margin_wall_chords. The domain size (Lx, Ly) is DERIVED from those
//   margins plus the grid footprint -- not the other way around. Every
//   airfoil's positioning is also bounded using a conservative rotation
//   radius (distance from pivot to the farther chord end), so the margin
//   guarantee holds regardless of what angle theta_ij ends up being --
//   this matters here because angles are meant to vary/be swept, unlike
//   the fixed-AoA single airfoil case.
// ============================================================================

#pragma once

#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <stdexcept>
#include <functional>
#include <algorithm>
#include <utility>

namespace geom {

typedef double T;
typedef long long plint;

// ---------------------------------------------------------------------
// NACA0012 half-thickness envelope, normalised chord x in [0,1]
// (identical formula to the one in NACA0012Palabos_Multigrid.cpp)
// ---------------------------------------------------------------------
inline T naca0012Thickness(T x) {
    if (x < 0.0) x = 0.0;
    if (x > 1.0) x = 1.0;
    return 5.0 * 0.12 * (
        0.2969 * std::sqrt(x)
        - 0.1260 * x
        - 0.3516 * x * x
        + 0.2843 * x * x * x
        - 0.1015 * x * x * x * x
    );
}

// ---------------------------------------------------------------------
// One airfoil unit in the lattice
// ---------------------------------------------------------------------
struct AirfoilUnit {
    int i = 0, j = 0;      // grid indices (row i, column j)
    T cx = 0.0, cy = 0.0;  // pivot location, COARSE lattice units
    T theta_rad = 0.0;     // rotation angle about the pivot, +ve = nose up
};

// ---------------------------------------------------------------------
// Lattice configuration. Distances that describe layout are given in
// CHORD units (not lattice cells) so the layout is resolution-independent
// -- change N_chord for a finer mesh and the geometry stays the same
// shape, exactly like x_foil/Lx scale with N_chord in the solid solver.
// ---------------------------------------------------------------------
struct LatticeParams {
    int nx = 6;               // airfoils per row (columns)
    int ny = 4;                // airfoils per column (rows)
    plint N_chord = 200;         // chord length, COARSE lattice units

    T pitch_x_chords = 2.5;     // horizontal pivot-to-pivot spacing, chords
    T pitch_y_chords = 2.5;     // vertical pivot-to-pivot spacing, chords

    T margin_inlet_chords  = 8.0;   // clearance ahead of column 0
    T margin_outlet_chords = 14.0;  // clearance behind the last column
    T margin_wall_chords   = 5.0;   // clearance above top row / below bottom row

    T pivot_chords = 0.25;      // rotation pivot on the chord (0=LE, 0.25=quarter-chord, 1=TE)

    // theta_deg_fn(i, j) supplies each airfoil's angle in degrees.
    // Default: uniform 0 deg (a flat, maximally transparent lattice).
    // Swap in a lookup table, a smooth field, or an optimiser's output
    // to get the "dynamic porous media" behaviour.
    std::function<T(int, int)> theta_deg_fn =
        [](int, int) { return 0.0; };
};

// ---------------------------------------------------------------------
// Result of building a lattice: unit list + the domain size it requires
// ---------------------------------------------------------------------
struct LatticeGeometry {
    std::vector<AirfoilUnit> units;
    plint Lx = 0, Ly = 0;  // derived domain size, COARSE lattice units
};

// ---------------------------------------------------------------------
// Build the lattice. Computes every airfoil's pivot location and the
// domain size that keeps every outer airfoil the requested margin away
// from its nearest boundary, for ANY rotation angle.
// ---------------------------------------------------------------------
inline LatticeGeometry buildLattice(const LatticeParams& p) {
    if (p.nx < 1 || p.ny < 1) {
        throw std::invalid_argument("LatticeParams: nx and ny must be >= 1");
    }
    if (p.pivot_chords < 0.0 || p.pivot_chords > 1.0) {
        throw std::invalid_argument("LatticeParams: pivot_chords must be in [0,1]");
    }

    LatticeGeometry g;
    const T C = static_cast<T>(p.N_chord);

    // Conservative bounding radius from the pivot to the farther chord
    // end, plus a small allowance for thickness (~6% chord). This is the
    // worst-case reach of the airfoil at ANY rotation angle, so margins
    // measured from the pivot using this radius stay valid however
    // theta_ij is set.
    const T R = std::max(p.pivot_chords, 1.0 - p.pivot_chords) + 0.06;

    auto colPivot_chords = [&](int j) {
        return p.margin_inlet_chords + R + j * p.pitch_x_chords;
    };
    auto rowPivot_chords = [&](int i) {
        return p.margin_wall_chords + R + i * p.pitch_y_chords;
    };

    T lastColPivot = colPivot_chords(p.nx - 1);
    T lastRowPivot = rowPivot_chords(p.ny - 1);

    T Lx_chords = lastColPivot + R + p.margin_outlet_chords;
    T Ly_chords = lastRowPivot + R + p.margin_wall_chords;

    g.Lx = static_cast<plint>(std::round(Lx_chords * C));
    g.Ly = static_cast<plint>(std::round(Ly_chords * C));

    g.units.reserve(static_cast<size_t>(p.nx) * static_cast<size_t>(p.ny));
    for (int i = 0; i < p.ny; ++i) {
        T y_chords = rowPivot_chords(i);
        for (int j = 0; j < p.nx; ++j) {
            T x_chords = colPivot_chords(j);

            AirfoilUnit u;
            u.i = i;
            u.j = j;
            u.cx = x_chords * C;
            u.cy = y_chords * C;
            u.theta_rad = p.theta_deg_fn(i, j) * M_PI / 180.0;
            g.units.push_back(u);
        }
    }
    return g;
}

// ---------------------------------------------------------------------
// Point-in-airfoil test for ONE unit, in COARSE lattice coordinates.
// Rotates the query point into the airfoil's local, unrotated frame
// about its pivot, re-expresses it in chord-normalised LE coordinates,
// then applies the standard thickness test -- same test as
// isInsideAirfoilReal() in the solid solver, generalised to an arbitrary
// pivot and rotation.
// ---------------------------------------------------------------------
inline bool isInsideUnit(T x, T y, const AirfoilUnit& u, plint N_chord, T pivot_chords = 0.25) {
    T dx = x - u.cx;
    T dy = y - u.cy;

    // Rotate by -theta to go from world frame into the airfoil's own frame.
    T c = std::cos(-u.theta_rad), s = std::sin(-u.theta_rad);
    T dx_local = dx * c - dy * s;
    T dy_local = dx * s + dy * c;

    const T C = static_cast<T>(N_chord);
    T x_norm = dx_local / C + pivot_chords;
    T y_norm = dy_local / C;

    if (x_norm < 0.0 || x_norm > 1.0) return false;
    T yt = naca0012Thickness(x_norm);
    return std::fabs(y_norm) < yt;
}

// ---------------------------------------------------------------------
// Point-in-lattice test: true if (x,y) lies inside ANY airfoil in the
// lattice. A cheap bounding-circle reject skips the exact test for the
// large majority of cells that are nowhere near a given airfoil.
// ---------------------------------------------------------------------
inline bool isInsideLattice(T x, T y, const std::vector<AirfoilUnit>& units,
                             plint N_chord, T pivot_chords = 0.25) {
    const T C = static_cast<T>(N_chord);
    const T R = std::max(pivot_chords, 1.0 - pivot_chords) + 0.06;
    const T rejectRadius = R * C;
    for (const auto& u : units) {
        T dx = x - u.cx, dy = y - u.cy;
        if (dx * dx + dy * dy > rejectRadius * rejectRadius) continue;
        if (isInsideUnit(x, y, u, N_chord, pivot_chords)) return true;
    }
    return false;
}

// ---------------------------------------------------------------------
// I/O -- write lattice centres/angles (for plotting, QA, or for the
// solver to read back in), and each unit's surface polyline (upper +
// lower), matching the writeAirfoilGeometry() CSV convention already
// used in NACA0012Palabos_Multigrid.cpp.
// ---------------------------------------------------------------------
inline void writeLatticeCenters(const LatticeGeometry& g, const std::string& filename) {
    std::ofstream file(filename);
    file << "i,j,cx,cy,theta_deg\n";
    for (const auto& u : g.units) {
        file << u.i << "," << u.j << "," << u.cx << "," << u.cy << ","
             << (u.theta_rad * 180.0 / M_PI) << "\n";
    }
}

inline void writeLatticeSurfaces(const LatticeGeometry& g, plint N_chord,
                                  const std::string& filename, int nPoints = 100,
                                  T pivot_chords = 0.25) {
    std::ofstream file(filename);
    file << "i,j,x_upper,y_upper,x_lower,y_lower\n";
    const T C = static_cast<T>(N_chord);
    for (const auto& u : g.units) {
        T c = std::cos(u.theta_rad), s = std::sin(u.theta_rad);
        for (int k = 0; k < nPoints; ++k) {
            T theta = M_PI * static_cast<T>(k) / static_cast<T>(nPoints - 1);
            T x_norm = 0.5 * (1.0 - std::cos(theta));
            T yt = naca0012Thickness(x_norm);

            // Local (unrotated) coords relative to the pivot.
            T xl = (x_norm - pivot_chords) * C;
            T yl_u = yt * C;
            T yl_l = -yt * C;

            T xu = u.cx + xl * c - yl_u * s;
            T yu = u.cy + xl * s + yl_u * c;
            T xd = u.cx + xl * c - yl_l * s;
            T yd = u.cy + xl * s + yl_l * c;

            file << u.i << "," << u.j << "," << xu << "," << yu << ","
                 << xd << "," << yd << "\n";
        }
    }
}

// ---------------------------------------------------------------------
// "Treat the whole lattice as one object" margin sizing.
//
// Standard bluff-body CFD practice scales inlet/outlet/wall clearance
// against the BODY's own size, not an arbitrary length. For a single
// airfoil, chord IS the body size, so margin_*_chords was already doing
// this. For a lattice, the relevant body is the WHOLE ARRAY -- so this
// computes the array's own streamwise/spanwise footprint and scales
// margins off THAT instead of a single chord.
//
// Usage: set nx, ny, pitch_x_chords, pitch_y_chords, pivot_chords on your
// LatticeParams as normal, then call applyBodyScaledMargins(p, ...) to
// OVERWRITE margin_inlet_chords / margin_outlet_chords / margin_wall_chords
// with values proportional to the array's footprint. Call this BEFORE
// buildLattice().
// ---------------------------------------------------------------------
struct ArrayFootprint {
    T width_chords;   // full streamwise extent: LE of column 0 to TE of last column
    T height_chords;  // full spanwise extent: bottom of row 0 to top of last row
};

inline ArrayFootprint computeArrayFootprint(int nx, int ny, T pitch_x_chords, T pitch_y_chords,
                                             T pivot_chords = 0.25) {
    T R = std::max(pivot_chords, 1.0 - pivot_chords) + 0.06;
    ArrayFootprint f;
    f.width_chords  = (nx - 1) * pitch_x_chords + 2.0 * R;
    f.height_chords = (ny - 1) * pitch_y_chords + 2.0 * R;
    return f;
}

// inlet_widths / outlet_widths are multiples of the array's own WIDTH;
// wall_heights is a multiple of the array's own HEIGHT. Typical bluff-body
// starting points: inlet ~1.5x, outlet ~3x (wake needs more room), wall
// ~1.5x -- tune per your blockage-ratio tolerance.
inline void applyBodyScaledMargins(LatticeParams& p,
                                    T inlet_widths = 1.5,
                                    T outlet_widths = 3.0,
                                    T wall_heights = 1.5) {
    ArrayFootprint f = computeArrayFootprint(p.nx, p.ny, p.pitch_x_chords, p.pitch_y_chords, p.pivot_chords);
    p.margin_inlet_chords  = inlet_widths  * f.width_chords;
    p.margin_outlet_chords = outlet_widths * f.width_chords;
    p.margin_wall_chords   = wall_heights  * f.height_chords;
}

// ---------------------------------------------------------------------
// Overlap check -- as pitch tightens, this verifies (by sampling each
// airfoil's actual rotated surface, not just conservative bounding
// circles) that no two units actually intersect. Useful specifically
// BECAUSE tightening pitch trades away the conservative safety margin.
// ---------------------------------------------------------------------
inline bool unitsOverlap(const AirfoilUnit& a, const AirfoilUnit& b,
                          plint N_chord, T pivot_chords = 0.25, int nSamples = 60) {
    const T C = static_cast<T>(N_chord);
    const T R = std::max(pivot_chords, 1.0 - pivot_chords) + 0.06;
    T dx = a.cx - b.cx, dy = a.cy - b.cy;
    if (dx*dx + dy*dy > (2*R*C) * (2*R*C)) return false; // cheap reject

    for (int k = 0; k < nSamples; ++k) {
        T theta = M_PI * static_cast<T>(k) / static_cast<T>(nSamples - 1);
        T x_norm = 0.5 * (1.0 - std::cos(theta));
        T yt = naca0012Thickness(x_norm);
        T xl = (x_norm - pivot_chords) * C;
        for (int side = 0; side < 2; ++side) {
            T yl = (side == 0) ? yt * C : -yt * C;
            T ca = std::cos(a.theta_rad), sa = std::sin(a.theta_rad);
            T xw = a.cx + xl * ca - yl * sa;
            T yw = a.cy + xl * sa + yl * ca;
            if (isInsideUnit(xw, yw, b, N_chord, pivot_chords)) return true;
        }
    }
    return false;
}

inline std::vector<std::pair<int,int>> findOverlaps(const LatticeGeometry& g, plint N_chord, T pivot_chords = 0.25) {
    std::vector<std::pair<int,int>> overlaps;
    for (size_t a = 0; a < g.units.size(); ++a) {
        for (size_t b = a + 1; b < g.units.size(); ++b) {
            if (unitsOverlap(g.units[a], g.units[b], N_chord, pivot_chords)) {
                overlaps.push_back(std::make_pair(static_cast<int>(a), static_cast<int>(b)));
            }
        }
    }
    return overlaps;
}

} // namespace geom

// ============================================================================
// SOLVER HOOK (reference only -- not compiled here, no Palabos include in
// this file). This is the shape of the stamping class you'd add to
// NACA0012Palabos_Multigrid.cpp, mirroring SetAirfoilSolidOnLevel but
// backed by geom::isInsideLattice() instead of isInsideAirfoilReal():
//
//   class SetLatticeSolidOnLevel : public BoxProcessingFunctional2D_S<int> {
//   public:
//       SetLatticeSolidOnLevel(plint scale_, const geom::LatticeGeometry& g_,
//                               plint N_chord_)
//           : scale(scale_), g(g_), N_chord(N_chord_) { }
//       virtual void process(Box2D domain, ScalarField2D<int>& flags) {
//           Dot2D offset = flags.getLocation();
//           for (plint ix = domain.x0; ix <= domain.x1; ++ix) {
//               for (plint iy = domain.y0; iy <= domain.y1; ++iy) {
//                   plint gx = ix + offset.x, gy = iy + offset.y;
//                   T x_coarse = static_cast<T>(gx) / scale;
//                   T y_coarse = static_cast<T>(gy) / scale;
//                   if (geom::isInsideLattice(x_coarse, y_coarse, g.units, N_chord)) {
//                       flags.get(ix, iy) = cellFlag::SOLID;
//                   }
//               }
//           }
//       }
//       // clone() / getTypeOfModification() as in SetAirfoilSolidOnLevel
//   private:
//       plint scale;
//       geom::LatticeGeometry g;
//       plint N_chord;
//   };
//
// param::Lx / param::Ly in the solver would then be set FROM
// geom::buildLattice(...).Lx / .Ly instead of the fixed 36*N_chord /
// 16*N_chord single-airfoil formula, and forceIds would become a
// std::vector<Array<plint,2>>, one pair per unit, if you want per-airfoil
// loads rather than a single lattice-wide force.
// ============================================================================