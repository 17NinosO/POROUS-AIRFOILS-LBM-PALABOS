// ============================================================================
// NACA0012Palabos_Multigrid.cpp
//
// Multigrid (2-level) variant of NACA0012Palabos.cpp.
//
// DESIGN NOTE (read before editing numLevel):
//   param::N_chord is the COARSE/level-0 chord length in lattice units, same
//   role it always had in the single-level file. Everything downstream
//   (Lx, Ly, x_foil, y_foil, geometry normalisation, tau/omega/nu_lb) is
//   built from it exactly as before.
//
//   The refined (finest) level automatically doubles resolution per level,
//   so the EFFECTIVE resolution at the airfoil surface is:
//       N_chord * 2^(numLevel - 1)
//   We fix that effective/target resolution at N_chord_fine = 1024 and
//   derive N_chord (coarse) from it. Change numLevel and N_chord updates
//   automatically to keep the fine-level resolution at 1024.
//
//   omega/tau is computed once, at the coarse (behaviorLevel) resolution,
//   and passed into makeBulkDynamics(). Palabos's MultiGridGenerator2D
//   rescales dynamics for finer levels internally — this matches the
//   pattern used successfully in the earlier N=200, 2-level MRT run on
//   this same project, not a fresh assumption.
// ============================================================================

#include "palabos2D.h"  //Declarations:class interfaces, types
#include "palabos2D.hh" //Definitions: template implementations

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cctype>
#include <string>
#include <algorithm>
#include <deque>

#include "complexDynamics/mrtDynamics.h"
#include "boundaryCondition/bounceBackModels.h"
#include "complexDynamics/mrtDynamics.hh"

#include "AirfoilLatticeGeometry.h"

using namespace plb;    //Palabos namespace
using namespace std;

typedef double T;       //All floating point uses this alias
#define DESCRIPTOR descriptors::MRTD2Q9Descriptor  //2D, 9-velocity lattice

// File-scope storage for the conforming lattice, built ONCE in main()
// before any stamping happens. isInsideAirfoilReal() reads from these.
// Not the cleanest design, but the minimal-diff, lowest-risk way to swap
// the geometry predicate without touching the proven stamping/BC/force/
// convergence machinery below, under a tight deadline.
// MOVED here (after using namespace plb;) -- plint is plb::plint, and the
// original placement before that using-directive left it out of scope,
// causing the "'plint' does not name a type" build error.
namespace conforming {
    geom::LatticeGeometry g;
    geom::SpatialHashGrid grid;
    plint miniChordCoarse = 1;
    plint N_chord_big = 1;   // the conforming lattice's own chord, coarse units -- needed for Cd/Cl normalization
    bool ready = false;
}

//===================================================================================================
//Simulation Parameters
//===================================================================================================

namespace param {

    //Physical Parameters
    const T Re = 1000.0; //Reynolds number
    const T AoA_deg = 11.0; //Angle of attack [degrees] (default; overridable via argv[1])
    const T AoA_rad = AoA_deg * M_PI / 180.0; //Angle of attack [radians]

    //---------------------------------------------------------------------
    // MULTIGRID RESOLUTION SETUP
    //---------------------------------------------------------------------
    // N_chord_fine: desired EFFECTIVE chord resolution at the airfoil surface
    // (the number we actually validate against the paper's N=512/1024 cases).
    const plint N_chord_fine = 400;

    // numLevel: number of grid levels. 2 = one coarse + one refined level.
    // Start at 2 and validate before going further (matches the earlier
    // conservative-first approach used throughout this project).
    const plint numLevel = 2;

    // N_chord: COARSE (level-0) chord length, derived so that the finest
    // level lands exactly on N_chord_fine. This is the constant used by
    // every downstream formula (domain size, geometry normalisation,
    // tau/omega), exactly as in the single-level file.
    const plint N_chord = N_chord_fine / (1 << (numLevel - 1));

    //Lattice velocity
    //Must satisfy Ma_lb = U_lb / c_s << 1 for incompressible flow.
    const T U_lb = 0.05;

    //Derived LBM Parameters (computed at the COARSE/behaviorLevel resolution)
    const T cs_sq = 1.0 / 3.0;
    const T nu_lb = U_lb * N_chord / Re; //Kinematic viscosity in lattice units.
    const T tau = 0.5 + nu_lb / cs_sq;   //Relaxation time for BGK collision operator.
    const T omega = 1.0 / tau;           //Relaxation rate for BGK collision operator.

    //MRT relaxation rates (D2Q9) — NOT currently wired into MRTdynamics'
    //constructor call below (which only takes omega). Kept here as
    //documented defaults for a future custom-rate MRT implementation;
    //don't assume these are actually affecting the simulation yet.
    const T s_1 = 1.4;
    const T s_2 = 1.4;
    const T s_4 = 1.2;
    const T s_6 = 1.2;

    //Domain Size [Lattice Units, COARSE level]
    //Same multiplier convention as the single-level file: 8C upstream,
    //28C downstream (36C total), 16C tall.
    const plint Lx = 36 * N_chord; // Domain Width  (coarse cells)
    const plint Ly = 16 * N_chord; // Domain Height (coarse cells)
    const plint x_foil = 8 * N_chord; // LE x position from the inlet.
    const plint y_foil = Ly / 2;      // LE y centred vertically.

    //Simulation Control Parameters
    const plint maxIter = 100000;
    const plint outIter = 1000;
    const T convTol = 1e-6; // relative fluctuation (std/|mean|) threshold for Cl AND Cd

    // Cp gets its OWN, LOOSER tolerance -- deliberately, not an oversight.
    // Cl/Cd are integrated quantities (summed momentum exchange over the
    // whole surface), which smooths out a lot of local noise for free.
    // Cp at a single station is a single lattice cell's density sample --
    // inherently noisier. Using the tight force tolerance here would
    // likely never be satisfied even once the flow is genuinely settled,
    // defeating the point of an early-stop check.
    const T cpConvTol = 1e-4;
    const plint forceLogIter = 10;

    // Surface Cp and the convergence check both use genuine ROLLING
    // windows (a circular buffer of the most recent samples), not "only
    // sample near a fixed endpoint" -- so whenever the run actually stops
    // (early convergence OR hitting maxIter), a valid trailing window is
    // always ready. cpAvgWindow doubles as the convergence-check window,
    // to avoid two near-duplicate "how much recent history" parameters.
    const plint cpAvgWindow = 5000;      // rolling window length, iterations
    const plint cpSampleIter = 100;      // Cp sampling cadence within the window
                                          // (coarser than forceLogIter -- each
                                          // sample costs 2*nStations MPI reductions)

    // Convergence check: every convCheckIter iterations (once past the
    // convMinIter warmup), evaluate relative fluctuation (std/|mean|) of
    // Cl and Cd over the last cpAvgWindow iterations. Must be under
    // convTol for convRequiredChecks CONSECUTIVE evaluations before
    // declaring convergence and stopping early -- a single quiet dip
    // isn't enough, avoiding the false-positive failure mode of the
    // earlier domain-averaged-energy check this replaces.
    const plint convCheckIter = 1000;
    const plint convMinIter = 20000;     // no convergence check before this many iterations
    const int convRequiredChecks = 3;

    // Collision model
    enum CollisionModel { BGK, MRT };
    const CollisionModel collisionModel = BGK;

    //Porosity Parameters (unused here, kept for interface parity with the
    //single-level file / porous case)
    const bool porous = false;
    const plint overlapWidth = 1;
    const plint behaviorLevel = 0;
    const T pore_width = 0.02;
    const int n_pores = 4;
    const T pore_centre_x = 0.30;
    const T pore_spacing = 0.08;

} //namespace param

//===================================================================================================
//Airfoil Geometry  (unchanged from single-level file — geometry is always
//expressed in level-0/coarse-equivalent coordinates, regardless of which
//actual grid level a cell being tested belongs to; see SetAirfoilSolidOnLevel)
//===================================================================================================

T naca0012Thickness(T x) {
    if (x<0.0) x=0.0;
    if (x>1.0) x=1.0;
    return 5.0 * 0.12 * (
        0.2969 * std::sqrt(x)
        - 0.1260 * x
        -0.3516 * x*x
        + 0.2843 * x*x*x
        - 0.1015 * x*x*x*x
    );
}

bool isInsidePore(plint ix, plint iy) {
    using namespace param;
    if (!porous) return false;
    T x_norm = static_cast<T>(ix - x_foil) / static_cast<T>(N_chord);
    for (int p = 0; p< n_pores; ++p) {
        T offset    =(static_cast<T>(p) - (n_pores - 1) / 2.0) * pore_spacing;
        T pore_x    = pore_centre_x + offset;
        if (std::fabs(x_norm - pore_x) < pore_width / 2.0){
            return true;
        }
    }
    return false;
}

// SWAPPED: was a solid-thickness-envelope test; now tests against the
// conforming mini-airfoil lattice (built once, before stamping, in
// initConformingGeometry() -- see main()). x_foil/y_foil/N_chord keep
// their EXACT original meaning and values -- the conforming lattice's own
// coordinates are LE-relative/vertically-centered, so translating by the
// same x_foil/y_foil offset the solid case always used places it in
// exactly the same domain location a solid airfoil of this chord would
// occupy. isInsidePore is now unreachable (porous=false by default) but
// left intact rather than removed, to minimize the diff under deadline
// pressure.
bool isInsideAirfoilReal(T x_coarse, T y_coarse) {
    using namespace param;
    T x_local = x_coarse - x_foil;
    T y_local = y_coarse - y_foil;
    std::vector<int> scratch; // thread-safe: local to this call, no shared state
    return geom::isInsideLatticeFast(x_local, y_local, conforming::g.units,
                                      conforming::grid, conforming::miniChordCoarse, scratch, 0.25);
}

// Builds the conforming lattice ONCE, before any level gets stamped.
// MUST be called before the stampAirfoilOnLevel loop in main().
void initConformingGeometry() {
    using namespace param;
    geom::ConformingLatticeParams cp;
    cp.N_chord_big = N_chord * 10;       // 10x bigger -- domain (Lx,Ly) UNCHANGED,
                                          // refine box now hugs actual geometry (fixed above)
                                          // instead of scaling with N_chord_big
    cp.c_mini_frac = 0.0246154;         // solved for ~30 columns at this N_chord_big -- SAME domain cost as 15/8 columns
                                         // SAME domain cost as 8 columns, but shows real
                                         // multi-row taper instead of degenerating to 1
                                         // row everywhere (verified visually before this change)
    cp.pitch_s_factor = 1.3;
    cp.pitch_n_factor = 1.0;
    cp.x_margin = 0.02;
    cp.pivot_chords = 0.25;

    conforming::g = geom::buildConformingLattice(cp);
    conforming::miniChordCoarse = geom::conformingChordUnits(cp);
    conforming::grid = geom::buildSpatialHash(conforming::g.units, conforming::miniChordCoarse, cp.pivot_chords);
    conforming::ready = true;
    conforming::N_chord_big = cp.N_chord_big;

    auto overlaps = geom::findOverlaps(conforming::g, conforming::miniChordCoarse, cp.pivot_chords);
    pcout << "Conforming lattice: " << conforming::g.units.size() << " mini-airfoils, "
          << "mini-chord=" << conforming::miniChordCoarse << " (coarse), "
          << "overlaps=" << overlaps.size() << "\n";
    if (!overlaps.empty()) {
        pcout << "WARNING: geometry has overlaps -- results will not be physical.\n";
    }

    geom::writeLatticeCenters(conforming::g, "conforming_centers.csv");
    geom::writeLatticeSurfaces(conforming::g, conforming::miniChordCoarse, "conforming_surfaces.csv");
}

namespace cellFlag {
    const int FLUID = 0;
    const int SOLID = 1;
    const int INLET = 2;
    const int OUTLET = 3;
    const int WALL = 4;
}

class SetAirfoilSolidOnLevel : public BoxProcessingFunctional2D_S<int> {
public:
    explicit SetAirfoilSolidOnLevel(plint scale_) : scale(scale_) { }

    virtual void process(Box2D domain, ScalarField2D<int>& flags) {
        Dot2D offset = flags.getLocation();
        for (plint ix = domain.x0; ix <= domain.x1; ++ix) {
            for (plint iy = domain.y0; iy <= domain.y1; ++iy) {
                plint gx = ix + offset.x;
                plint gy = iy + offset.y;
                T x_coarse = static_cast<T>(gx) / scale;
                T y_coarse = static_cast<T>(gy) / scale;
                if (isInsideAirfoilReal(x_coarse, y_coarse)) {
                    flags.get(ix, iy) = cellFlag::SOLID;
                }
            }
        }
    }

    virtual SetAirfoilSolidOnLevel* clone() const {
        return new SetAirfoilSolidOnLevel(*this);
    }

    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
    }

private:
    plint scale;
};

// CHANGED vs single-level file: takes a `trackForce` flag. Only the finest
// level gets MomentumExchangeBounceBack (force-accumulating); every other
// level gets plain BounceBack, so the airfoil footprint under the fine
// patch stays geometrically consistent (matches the earlier working
// multi-level runs on this project) WITHOUT double-counting force into
// forceIds from more than one level.
void stampAirfoilOnLevel(MultiGridLattice2D<T, DESCRIPTOR>& lattice, plint iLevel,
                          Array<plint,2> forceIds, bool trackForce) {
    MultiBlockLattice2D<T, DESCRIPTOR>& comp = lattice.getComponent(iLevel);
    Box2D bb = comp.getBoundingBox();
    plint scale = 1 << iLevel;

    MultiScalarField2D<int> flags(comp);
    applyProcessingFunctional(new SetAirfoilSolidOnLevel(scale), bb, flags);

    if (trackForce) {
        defineDynamics(comp, flags, new MomentumExchangeBounceBack<T, DESCRIPTOR>(forceIds), cellFlag::SOLID);
    } else {
        defineDynamics(comp, flags, new BounceBack<T, DESCRIPTOR>(), cellFlag::SOLID);
    }
    pcout << "Airfoil stamped on level " << iLevel
          << (trackForce ? " (force-tracked)" : " (plain, no force tracking)") << "\n";
}

//===================================================================================================
//Surface Point Generator (unchanged)
//===================================================================================================

void writeAirfoilGeometry(const std::string& filename, int nPoints=200) {
    using namespace param;
    std::ofstream file(filename);
    file << "x_upper, y_upper, x_lower, y_lower\n";
    for (int i = 0; i < nPoints; ++i) {
        T theta = M_PI * static_cast<T>(i) / static_cast<T>(nPoints - 1);
        T x_norm = 0.5 * (1.0 - std::cos(theta));
        T yt = naca0012Thickness(x_norm);
        T x_lat = x_foil + x_norm * N_chord;
        T y_upper = y_foil + yt * N_chord;
        T y_lower = y_foil - yt * N_chord;
        file << x_lat << "," << y_upper << ","
            << x_lat << "," << y_lower << "\n";
    }
}

//===================================================================================================
//Dynamics factory (unchanged) — omega is built at the coarse/behaviorLevel
//resolution; Palabos rescales it for finer levels internally.
//===================================================================================================

Dynamics<T, DESCRIPTOR>* makeBulkDynamics(T omega_runtime, param::CollisionModel model) {
    if (model == param::BGK)
        return new BGKdynamics<T, DESCRIPTOR>(omega_runtime);
    else
        return new MRTdynamics<T, DESCRIPTOR>(omega_runtime);
}

//===================================================================================================
//Initialisation (unchanged)
//===================================================================================================

void initializeDomain(MultiGridLattice2D<T, DESCRIPTOR>& lattice, T AoA_rad) {
    using namespace param;
    const T rho_0 = 1.0;
    Array<T,2> u_inf(U_lb * std::cos(AoA_rad), U_lb * std::sin(AoA_rad));

    for (plint iLevel = 0; iLevel < lattice.getNumLevels(); ++iLevel) {
        initializeAtEquilibrium(
            lattice.getComponent(iLevel),
            lattice.getComponent(iLevel).getBoundingBox(),
            rho_0, u_inf);
    }
    lattice.initialize();
}

//===================================================================================================
//Boundary Conditions (unchanged — still applied on the coarse/level-0 block
//only, which is correct: domain-edge boundaries always live on level 0)
//===================================================================================================

void setBoundaryConditions(MultiGridLattice2D<T, DESCRIPTOR>& lattice, T AoA_rad) {
    using namespace param;
    Array<T, 2> u_inf(U_lb * std::cos(AoA_rad), U_lb * std::sin(AoA_rad));

    MultiBlockLattice2D<T, DESCRIPTOR>& coarse = lattice.getComponent(0);

    OnLatticeBoundaryCondition2D<T, DESCRIPTOR>* bc =
        createLocalBoundaryCondition2D<T, DESCRIPTOR>();

    // Inlet (left) — velocity
    bc->addVelocityBoundary0N(Box2D(0, 0, 1, Ly-2), coarse);
    setBoundaryVelocity(coarse, Box2D(0, 0, 1, Ly-2), u_inf);

    // Outlet (right) — pressure
    bc->addPressureBoundary0P(Box2D(Lx-1, Lx-1, 1, Ly-2), coarse);
    setBoundaryDensity(coarse, Box2D(Lx-1, Lx-1, 1, Ly-2), (T)1.0);

    // Top and bottom — velocity, matching free-stream (Gabriel's reference approach)
    bc->addVelocityBoundary1P(Box2D(1, Lx-2, Ly-1, Ly-1), coarse);
    setBoundaryVelocity(coarse, Box2D(1, Lx-2, Ly-1, Ly-1), u_inf);
    bc->addVelocityBoundary1N(Box2D(1, Lx-2, 0, 0), coarse);
    setBoundaryVelocity(coarse, Box2D(1, Lx-2, 0, 0), u_inf);

    defineDynamics(coarse, Box2D(0,0,0,0),       new BounceBack<T,DESCRIPTOR>());
    defineDynamics(coarse, Box2D(0,0,Ly-1,Ly-1), new BounceBack<T,DESCRIPTOR>());
    defineDynamics(coarse, Box2D(Lx-1,Lx-1,0,0), new BounceBack<T,DESCRIPTOR>());
    defineDynamics(coarse, Box2D(Lx-1,Lx-1,Ly-1,Ly-1), new BounceBack<T,DESCRIPTOR>());

    delete bc;
}

//===================================================================================================
//Force Computation (unchanged)
//===================================================================================================

void computeCoefficients(T Fx, T Fy, T& Cl, T& Cd, T AoA_rad) {
    using namespace param;
    T dx = 1.0 / (T)N_chord_fine;       // force measured on finest level — use FINE chord length
    T dt = U_lb / (T)N_chord_fine;      // (assumes convective scaling: dx/dt ratio preserved across levels)
    T scale = dx*dx*dx / (dt*dt);

    T Fx_phys = 2.0 * Fx * scale;
    T Fy_phys = 2.0 * Fy * scale;

    T q_inf = 0.5 * 1.0 * 1.0;
    // FIXED: was area=1.0, a leftover from the original single-solid-airfoil
    // file where the airfoil's own chord WAS the "1" in this physical unit
    // system (dx=1/N_chord_fine sets 1 physical unit = 1 original chord).
    // The conforming lattice's own chord is now N_chord_big (coarse units),
    // which is N_chord_big/N_chord times bigger than that original
    // reference chord -- area must scale accordingly, or Cl/Cd silently
    // stay normalized against a chord that no longer describes this
    // geometry at all.
    T area = static_cast<T>(conforming::N_chord_big) / static_cast<T>(N_chord);

    T lift = -Fx_phys * std::sin(AoA_rad) + Fy_phys * std::cos(AoA_rad);
    T drag =  Fx_phys * std::cos(AoA_rad) + Fy_phys * std::sin(AoA_rad);

    Cl = lift / (q_inf * area);
    Cd = drag / (q_inf * area);
}

//===================================================================================================
//Surface Cp -- WINDOWED MEAN, not a single snapshot. A single-iteration
//sample only captures whichever phase of an oscillation (vortex
//shedding, unsteady separation) happened to land on that exact
//timestep -- especially relevant since this file has no convergence-
//based early stop, so every run goes the full maxIter regardless of
//whether the flow ever truly settles. This accumulates mean + std at
//each station over a trailing window of iterations, same philosophy as
//the windowed-mean force statistics already used elsewhere in this
//project for high-AoA unsteady cases.
//
// Samples at fine-level cells just OUTSIDE the solid surface (offset by
// a few cells along y -- the airfoil is axis-aligned in domain
// coordinates, AoA is applied via inflow direction not geometry
// rotation, so "outside the surface" is a vertical offset here).
// Converts density to Cp via p = rho * cs^2. Uses computeAverageDensity()
// on 1-cell boxes -- MPI-safe regardless of domain decomposition, unlike
// a raw lattice.get() call.
//===================================================================================================

class SurfaceCpAccumulator {
public:
    SurfaceCpAccumulator(int nStations_, size_t maxSamples_) : nStations(nStations_),
        maxSamples(maxSamples_), bufUpper(nStations_), bufLower(nStations_), consecutivePasses(0) { }

    // FIXED: was sampling along the OLD solid-airfoil's surface formula
    // (x_foil + x_norm*N_chord), which describes where a single solid
    // airfoil of the file's original fixed chord would sit -- completely
    // wrong location once the geometry became 87 separate, individually
    // rotated mini-airfoils spanning a much larger N_chord_big. Now
    // samples each mini-airfoil's own MID-CHORD point (upper + lower),
    // using its actual position/rotation from conforming::g.units, with
    // the sample offset applied in the airfoil's LOCAL frame (before
    // rotating) so "outward" is correct regardless of that airfoil's
    // individual angle -- not just a world-frame +/-y offset, which would
    // only be correct for an axis-aligned airfoil.
    void accumulate(MultiGridLattice2D<T, DESCRIPTOR>& lattice, plint interpLevel) {
        using namespace param;
        const plint scale = 1 << interpLevel;
        const T q_inf_lattice = 0.5 * U_lb * U_lb;
        const plint sampleOffsetCells = 3;
        const T C = static_cast<T>(conforming::miniChordCoarse);
        const T x_norm = 0.5;   // mid-chord: thickest, most representative station
        const T yt = naca0012Thickness(x_norm);
        const T xl = (x_norm - 0.25) * C;   // pivot = quarter-chord, matches geometry convention

        for (int k = 0; k < nStations; ++k) {
            const geom::AirfoilUnit& u = conforming::g.units[k];
            T c = std::cos(u.theta_rad), s = std::sin(u.theta_rad);

            // Offset outward in the LOCAL (unrotated) frame, THEN rotate --
            // this is what makes "outward" correct for each unit's own angle.
            T yl_u = yt * C + sampleOffsetCells;
            T yl_l = -yt * C - sampleOffsetCells;

            T xu_local = u.cx + xl*c - yl_u*s, yu_local = u.cy + xl*s + yl_u*c;
            T xd_local = u.cx + xl*c - yl_l*s, yd_local = u.cy + xl*s + yl_l*c;

            // Translate LE-relative lattice coords into world coords --
            // same offset isInsideAirfoilReal() uses.
            plint ix_upper_fine = static_cast<plint>(std::round((xu_local + x_foil) * scale));
            plint iy_upper_fine = static_cast<plint>(std::round((yu_local + y_foil) * scale));
            plint ix_lower_fine = static_cast<plint>(std::round((xd_local + x_foil) * scale));
            plint iy_lower_fine = static_cast<plint>(std::round((yd_local + y_foil) * scale));

            T rho_upper = computeAverageDensity(lattice.getComponent(interpLevel),
                                                 Box2D(ix_upper_fine, ix_upper_fine, iy_upper_fine, iy_upper_fine));
            T rho_lower = computeAverageDensity(lattice.getComponent(interpLevel),
                                                 Box2D(ix_lower_fine, ix_lower_fine, iy_lower_fine, iy_lower_fine));

            T Cp_upper = (rho_upper - 1.0) * cs_sq / q_inf_lattice;
            T Cp_lower = (rho_lower - 1.0) * cs_sq / q_inf_lattice;

            // Rolling: drop the oldest sample once at capacity, so this
            // always reflects the most recent window regardless of when
            // the run actually stops.
            bufUpper[k].push_back(Cp_upper);
            if (bufUpper[k].size() > maxSamples) bufUpper[k].pop_front();
            bufLower[k].push_back(Cp_lower);
            if (bufLower[k].size() > maxSamples) bufLower[k].pop_front();
        }
    }

    // Mirrors ForceConvergenceTracker::checkConverged()'s interface, so the
    // two checks combine symmetrically in the main loop. Uses the WORST-CASE
    // (max) relative fluctuation across every station and both surfaces --
    // not an average -- because a single still-noisy station (e.g. near a
    // separation point) is exactly the failure mode this check exists to
    // catch, and averaging would let many quiet stations hide it.
    bool checkConverged(T tol, int requiredConsecutive) {
        if (nStations == 0 || bufUpper[0].size() < maxSamples) {
            consecutivePasses = 0;
            pcout << "  [Cp convergence check] window filling: "
                  << (nStations > 0 ? bufUpper[0].size() : 0) << "/" << maxSamples << " samples\n";
            return false;
        }

        T maxRelFluct = 0.0;
        for (int k = 0; k < nStations; ++k) {
            T meanU, stdU, meanL, stdL;
            statsOf(bufUpper[k], meanU, stdU);
            statsOf(bufLower[k], meanL, stdL);
            T relU = std::fabs(meanU) > 1e-6 ? stdU / std::fabs(meanU) : stdU;
            T relL = std::fabs(meanL) > 1e-6 ? stdL / std::fabs(meanL) : stdL;
            maxRelFluct = std::max(maxRelFluct, std::max(relU, relL));
        }

        bool passed = maxRelFluct < tol;
        consecutivePasses = passed ? (consecutivePasses + 1) : 0;

        pcout << "  [Cp convergence check] max_relFluct=" << maxRelFluct
              << " (tol=" << tol << ") consecutive_passes=" << consecutivePasses
              << "/" << requiredConsecutive << "\n";

        return consecutivePasses >= requiredConsecutive;
    }

    void writeToFile(const std::string& filename) const {
        std::ofstream file(filename);
        // RELABELED: each station is now one MINI-AIRFOIL (unit_i, unit_j
        // from the conforming lattice), not a fractional x/c position
        // along one continuous surface -- the old "x_over_c" label no
        // longer means anything here.
        file << "unit_i,unit_j,x_world_coarse,Cp_upper_mean,Cp_upper_std,Cp_lower_mean,Cp_lower_std,n_samples\n";
        for (int k = 0; k < nStations; ++k) {
            const geom::AirfoilUnit& u = conforming::g.units[k];
            T x_world = u.cx + param::x_foil;
            T meanU, stdU, meanL, stdL;
            statsOf(bufUpper[k], meanU, stdU);
            statsOf(bufLower[k], meanL, stdL);
            file << u.i << "," << u.j << "," << x_world << "," << meanU << "," << stdU << ","
                 << meanL << "," << stdL << "," << bufUpper[k].size() << "\n";
        }
        pcout << "Windowed-mean surface Cp (" << (nStations > 0 ? bufUpper[0].size() : 0)
              << " samples) written to " << filename << "\n";
    }

private:
    static void statsOf(const std::deque<T>& buf, T& mean, T& stddev) {
        if (buf.empty()) { mean = 0.0; stddev = 0.0; return; }
        T sum = 0.0, sumSq = 0.0;
        for (T v : buf) { sum += v; sumSq += v * v; }
        mean = sum / buf.size();
        T var = sumSq / buf.size() - mean * mean;
        stddev = var > 0.0 ? std::sqrt(var) : 0.0;
    }

    int nStations;
    size_t maxSamples;
    std::vector<std::deque<T>> bufUpper, bufLower;
    int consecutivePasses;
};

//===================================================================================================
//Force-based convergence check -- replaces the old domain-averaged-energy
//ValueTracer (removed earlier for being prone to false positives: energy
//averaged over the whole domain can plateau even while the near-airfoil
//flow is still genuinely transient). Tracks Cl and Cd directly in a
//rolling window -- the quantities that actually matter for this
//comparison -- and requires the relative fluctuation (std/|mean|) of
//BOTH to stay under tolerance for several CONSECUTIVE checks before
//declaring convergence, not just one quiet moment.
//===================================================================================================

class ForceConvergenceTracker {
public:
    ForceConvergenceTracker(size_t maxSamples_) : maxSamples(maxSamples_), consecutivePasses(0) { }

    void addSample(T Cl, T Cd) {
        clBuf.push_back(Cl);
        if (clBuf.size() > maxSamples) clBuf.pop_front();
        cdBuf.push_back(Cd);
        if (cdBuf.size() > maxSamples) cdBuf.pop_front();
    }

    // Returns true only once convRequiredChecks consecutive calls have
    // all passed the tolerance test. Call this once per convCheckIter.
    bool checkConverged(T tol, int requiredConsecutive) {
        if (clBuf.size() < maxSamples) {
            consecutivePasses = 0;
            pcout << "  [convergence check] window filling: " << clBuf.size() << "/" << maxSamples << " samples\n";
            return false;
        }

        T clMean, clStd, cdMean, cdStd;
        statsOf(clBuf, clMean, clStd);
        statsOf(cdBuf, cdMean, cdStd);

        T clRelFluct = std::fabs(clMean) > 1e-12 ? clStd / std::fabs(clMean) : clStd;
        T cdRelFluct = std::fabs(cdMean) > 1e-12 ? cdStd / std::fabs(cdMean) : cdStd;

        bool passed = (clRelFluct < tol) && (cdRelFluct < tol);
        consecutivePasses = passed ? (consecutivePasses + 1) : 0;

        pcout << "  [convergence check] Cl_relFluct=" << clRelFluct
              << " Cd_relFluct=" << cdRelFluct << " (tol=" << tol << ") "
              << "consecutive_passes=" << consecutivePasses << "/" << requiredConsecutive << "\n";

        return consecutivePasses >= requiredConsecutive;
    }

private:
    static void statsOf(const std::deque<T>& buf, T& mean, T& stddev) {
        T sum = 0.0, sumSq = 0.0;
        for (T v : buf) { sum += v; sumSq += v * v; }
        mean = sum / buf.size();
        T var = sumSq / buf.size() - mean * mean;
        stddev = var > 0.0 ? std::sqrt(var) : 0.0;
    }

    size_t maxSamples;
    std::deque<T> clBuf, cdBuf;
    int consecutivePasses;
};

// Writes a 1/0 field marking which cells are solid airfoil material, so
// the actual geometry is directly visible in ParaView -- previously
// writeVTK() only wrote flow quantities (velocity/density), with no way
// to see WHERE the mini-airfoils are except by inferring it from flow
// disturbance, which is nearly invisible this early in a run. Reuses
// isInsideAirfoilReal() directly (same predicate the stamping functional
// uses), so this is guaranteed to match the actual stamped geometry.
class ComputeSolidMask2D : public BoxProcessingFunctional2D_S<T> {
public:
    explicit ComputeSolidMask2D(plint scale_) : scale(scale_) { }

    virtual void process(Box2D domain, ScalarField2D<T>& mask) {
        Dot2D offset = mask.getLocation();
        for (plint ix = domain.x0; ix <= domain.x1; ++ix) {
            for (plint iy = domain.y0; iy <= domain.y1; ++iy) {
                plint gx = ix + offset.x, gy = iy + offset.y;
                T x_coarse = static_cast<T>(gx) / scale;
                T y_coarse = static_cast<T>(gy) / scale;
                mask.get(ix, iy) = isInsideAirfoilReal(x_coarse, y_coarse) ? 1.0 : 0.0;
            }
        }
    }

    virtual ComputeSolidMask2D* clone() const { return new ComputeSolidMask2D(*this); }

    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
    }

private:
    plint scale;
};

// CHANGED vs single-level file: level-0 write is now a much bigger physical
// domain at coarser spacing (still full Lx x Ly, coarse cells). Level>0
// cropping to refineBox*scale is unchanged. ADDED: isSolid field so the
// mini-airfoil geometry is directly visible in ParaView, not just inferred.
void writeVTK(MultiGridLattice2D<T, DESCRIPTOR>& lattice, plint iter, Box2D const& refineBox) {
    for (plint iLevel = 0; iLevel < lattice.getNumLevels(); ++iLevel) {
        MultiBlockLattice2D<T, DESCRIPTOR>& c = lattice.getComponent(iLevel);

        Box2D writeBox = c.getBoundingBox();
        if (iLevel > 0) {
            plint scale = 1 << iLevel;
            writeBox = Box2D(refineBox.x0 * scale, refineBox.x1 * scale,
                              refineBox.y0 * scale, refineBox.y1 * scale);
        }

        MultiScalarField2D<T> solidMask(c);
        applyProcessingFunctional(new ComputeSolidMask2D(1 << iLevel), writeBox, solidMask);

        VtkImageOutput2D<T> vtkOut(
            createFileName("vtk_level" + std::to_string(iLevel), iter, 6),
            1.0 / std::pow(2.0, iLevel));
        vtkOut.writeData<float>(*computeVelocityNorm(c, writeBox), "velocityNorm", 1.0);
        vtkOut.writeData<2, float>(*computeVelocity(c, writeBox), "velocity", 1.0);
        vtkOut.writeData<float>(*computeDensity(c, writeBox), "density", 1.0);
        vtkOut.writeData<float>(solidMask, "isSolid", 1.0);
    }
}

//===================================================================================================
//Simulation Loop (unchanged from the flush/header-fixed single-level file)
//===================================================================================================
void runSimulation(MultiGridLattice2D<T, DESCRIPTOR>& lattice, Box2D const& refineBox,
                    plint interpLevel, Array<plint,2> forceIds, T AoA_rad, const std::string& runTag)
{
    using namespace param;
    std::string forcesFilename = "forces_" + runTag + ".txt";
    std::ofstream forceFile(forcesFilename);
    forceFile << "iter,Cl,Cd\n";

    // Both rolling windows sample continuously from iteration 0 -- NOT
    // "only near a fixed endpoint" -- so whichever way this run ends
    // (early convergence or hitting maxIter), a valid trailing window is
    // always ready. cpAvgWindow/cpSampleIter and cpAvgWindow/forceLogIter
    // set each buffer's capacity in SAMPLES (not iterations).
    // Station count now matches the ACTUAL mini-airfoil count (one
    // mid-chord sample pair per unit), not a fixed 100 -- was previously
    // hardcoded independent of how many units the geometry actually has.
    SurfaceCpAccumulator cpAccum(static_cast<int>(conforming::g.units.size()),
                                  static_cast<size_t>(cpAvgWindow / cpSampleIter));
    ForceConvergenceTracker convTracker(static_cast<size_t>(cpAvgWindow / forceLogIter));

    pcout << "Starting simulation: " << maxIter << " max iterations (hard cap)\n";
    pcout << "Early stop requires BOTH checks to pass, " << convRequiredChecks
          << " consecutive times, no earlier than iter " << convMinIter << ":\n";
    pcout << "  1) Force check:  Cl/Cd relative fluctuation < " << convTol
          << " over a rolling " << cpAvgWindow << "-iteration window\n";
    pcout << "  2) Cp check:     max relative fluctuation across all surface stations < " << cpConvTol
          << " over the same window (sampled every " << cpSampleIter << " iterations)\n";

    plint finalIter = maxIter;
    bool convergedEarly = false;

    for (plint iT = 0; iT <= maxIter; ++iT) {

        if (iT % forceLogIter == 0) {
            T Fx = lattice.getComponent(interpLevel).getInternalStatistics().getSum(forceIds[0]);
            T Fy = lattice.getComponent(interpLevel).getInternalStatistics().getSum(forceIds[1]);

            // The fine level's collideAndStream() runs 2^interpLevel times per single
            // outer call (Palabos multigrid convective sub-stepping), but evaluateStatistics()
            // only fires once at the end — so getSum() returns the accumulated total across
            // all sub-steps, not one. Divide back down to a per-coarse-iteration value.
            const T subStepsPerCoarseIter = static_cast<T>(1 << interpLevel);
            Fx /= subStepsPerCoarseIter;
            Fy /= subStepsPerCoarseIter;

            T Cl, Cd;
            computeCoefficients(Fx, Fy, Cl, Cd, AoA_rad);
            forceFile << iT << "," << Cl << "," << Cd << "\n";
            forceFile.flush();
            convTracker.addSample(Cl, Cd);
        }

        if (iT % outIter == 0) {
            T avgEnergy = computeAverageEnergy(lattice.getComponent(0));
            pcout << "iter=" << iT << " E=" << avgEnergy << "\n";
        }

        if (iT >= convMinIter && iT % cpSampleIter == 0) {
            cpAccum.accumulate(lattice, interpLevel);
        }

        if (iT >= convMinIter && iT % convCheckIter == 0) {
            // Evaluate BOTH checks every time, unconditionally -- if this used
            // `forceOK && cpAccum.checkConverged(...)` directly, short-circuit
            // evaluation would skip the Cp check whenever forceOK is false,
            // leaving its consecutive-pass streak stale instead of correctly
            // updating (and likely resetting) every cycle.
            bool forceOK = convTracker.checkConverged(convTol, convRequiredChecks);
            bool cpOK = cpAccum.checkConverged(cpConvTol, convRequiredChecks);

            if (forceOK && cpOK) {
                pcout << "Converged (forces AND Cp both stable) at iter=" << iT
                      << " -- stopping early.\n";
                finalIter = iT;
                convergedEarly = true;
                writeVTK(lattice, iT, refineBox);
                break;
            }
        }

        if (iT == maxIter) {
            writeVTK(lattice, iT, refineBox);
        }

        lattice.collideAndStream();
    }

    forceFile.close();
    pcout << "Simulation complete at iter=" << finalIter
          << (convergedEarly ? " (converged early)" : " (hit maxIter cap, did not converge)")
          << ". Forces saved to " << forcesFilename << "\n";

    cpAccum.writeToFile("surface_cp_" + runTag + ".csv");
}

//===================================================================================================
//Main
//===================================================================================================

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);

    // Force every std::cout write to flush IMMEDIATELY, not batch into
    // large blocks -- default C++ stdio behavior block-buffers output
    // whenever it's not connected to an interactive terminal (i.e.
    // whenever it's piped through tee, redirected to a file, or run
    // under PBS with -k oe). Without this, everything the program has
    // genuinely already printed -- including every convergence check --
    // can sit invisible in a buffer for a long time before appearing in
    // the file you're tailing, even though nothing is actually wrong.
    std::cout << std::unitbuf;

    T AoA_deg_runtime = (argc > 1) ? std::atof(argv[1]) : param::AoA_deg;
    T AoA_rad_runtime = AoA_deg_runtime * M_PI / 180.0;

    T Re_runtime = (argc > 2) ? std::atof(argv[2]) : param::Re;

    param::CollisionModel collisionModel_runtime = param::collisionModel;
    std::string modelStr = (param::collisionModel == param::BGK) ? "BGK" : "MRT";
    if (argc > 3) {
        std::string arg3 = argv[3];
        for (auto& c : arg3) c = std::toupper(c);
        if (arg3 == "MRT") { collisionModel_runtime = param::MRT; modelStr = "MRT"; }
        else if (arg3 == "BGK") { collisionModel_runtime = param::BGK; modelStr = "BGK"; }
        else {
            pcout << "WARNING: unrecognised collision model '" << argv[3]
                  << "', expected BGK or MRT -- falling back to default.\n";
        }
    }

    // Runtime LBM parameters, derived from Re_runtime -- NOT param::nu_lb /
    // param::tau / param::omega, which are only compile-time DEFAULTS now.
    T nu_lb_runtime = param::U_lb * param::N_chord / Re_runtime;
    T tau_runtime = 0.5 + nu_lb_runtime / param::cs_sq;
    T omega_runtime = 1.0 / tau_runtime;

    // Run tag for unique output filenames across a Re x AoA x model sweep
    // sharing a directory (forces_<tag>.txt, surface_cp_<tag>.csv).
    std::ostringstream tagStream;
    tagStream << "Conforming30col_Re" << static_cast<plint>(Re_runtime)
              << "_AoA" << AoA_deg_runtime << "_" << modelStr;
    std::string runTag = tagStream.str();

    pcout << "============================================\n";
    pcout << "   NACA 0012 LBM Simulation - MULTIGRID    \n";
    pcout << "   (MRT vs BGK comparison build)           \n";
    pcout << "============================================\n";
    pcout << "Re          = " << Re_runtime      << "\n";
    pcout << "AoA         = " << AoA_deg_runtime  << " deg\n";
    pcout << "numLevel    = " << param::numLevel  << "\n";
    pcout << "N_chord(coarse) = " << param::N_chord << " cells\n";
    pcout << "N_chord(fine, effective) = " << param::N_chord_fine << " cells\n";
    pcout << "Domain(coarse)  = " << param::Lx << " x " << param::Ly << " cells\n";
    pcout << "tau (coarse)    = " << tau_runtime      << "\n";
    pcout << "omega (coarse)  = " << omega_runtime    << "\n";
    pcout << "U_lb        = " << param::U_lb     << "\n";
    pcout << "Max iters   = " << param::maxIter  << "\n";
    pcout << "Collision   = " << modelStr << "\n";
    pcout << "Run tag     = " << runTag << "\n";
    pcout << "============================================\n";

    writeAirfoilGeometry("airfoil_geometry.csv");
    pcout << "Airfoil geometry written to airfoil_geometry.csv (reference solid outline, unused for stamping)\n";

    initConformingGeometry();   // MUST happen before the stamping loop below

    plint numLevel = param::numLevel, overlapWidth = param::overlapWidth, behaviorLevel = param::behaviorLevel;
    MultiGridManagement2D management(param::Lx, param::Ly, numLevel, overlapWidth, behaviorLevel);

    // FIXED: refineBox now computed from the CONFORMING LATTICE's actual
    // bounding box (conforming::g.units), not the old "+-4*N_chord"
    // formula inherited from the single-solid-airfoil file. That formula
    // assumed N_chord_big WAS the domain-sizing chord and needed generous
    // wake margin around a solid body -- growing N_chord_big under that
    // formula would blow the refine box past the domain edges. The
    // lattice's real geometric footprint is much smaller than that
    // formula assumes (thin, mostly-open structure), so this instead
    // hugs the ACTUAL extent with a fixed buffer -- letting N_chord_big
    // scale independently of domain size, same principle already used
    // for the rectangular lattice's refine box earlier in this project.
    T R_mini = (std::max(0.25, 1.0 - 0.25) + 0.06) * static_cast<T>(conforming::miniChordCoarse);
    T ext_xmin = 1e18, ext_xmax = -1e18, ext_ymin = 1e18, ext_ymax = -1e18;
    for (auto& u : conforming::g.units) {
        ext_xmin = std::min(ext_xmin, u.cx - R_mini); ext_xmax = std::max(ext_xmax, u.cx + R_mini);
        ext_ymin = std::min(ext_ymin, u.cy - R_mini); ext_ymax = std::max(ext_ymax, u.cy + R_mini);
    }
    const T padNear = 5.0 * conforming::miniChordCoarse;   // tight fixed buffer, all sides
    const T padWake = 20.0 * conforming::miniChordCoarse;  // extra downstream, for near-wake capture
    plint rb_x0 = static_cast<plint>(param::x_foil + ext_xmin - padNear);
    plint rb_x1 = static_cast<plint>(param::x_foil + ext_xmax + padNear + padWake);
    plint rb_y0 = static_cast<plint>(param::y_foil + ext_ymin - padNear);
    plint rb_y1 = static_cast<plint>(param::y_foil + ext_ymax + padNear);
    Box2D refineBox(std::max(rb_x0, (plint)1), std::min(rb_x1, param::Lx - 2),
                     std::max(rb_y0, (plint)1), std::min(rb_y1, param::Ly - 2));
    pcout << "Refine box (from actual geometry, not old formula): x[" << refineBox.x0 << "," << refineBox.x1
          << "] y[" << refineBox.y0 << "," << refineBox.y1 << "]\n";
    management.refine(0, refineBox);

    plint xTiles = global::mpi().getSize(), yTiles = 1, interpLevel = numLevel - 1;
    ParallellizeBySquares2D* parallelizer = new ParallellizeBySquares2D(
        management.getBulks(), management.getBoundingBox(interpLevel), xTiles, yTiles);
    management.parallelize(parallelizer);

    MultiGridLattice2D<T, DESCRIPTOR> lattice =
        *MultiGridGenerator2D<T, DESCRIPTOR>::createRefinedLatticeCubicInterpolationFiltering(
            management, makeBulkDynamics(omega_runtime, collisionModel_runtime), behaviorLevel);
    pcout << "Lattice built (numLevel=" << numLevel << ")\n";

    Box2D bb0 = lattice.getComponent(0).getBoundingBox();
    pcout << "Level 0 bbox: x[" << bb0.x0 << "," << bb0.x1 << "] y[" << bb0.y0 << "," << bb0.y1 << "]\n";
    
    Box2D bb1 = lattice.getComponent(interpLevel).getBoundingBox();
    pcout << "Level " << interpLevel << " bbox: x[" << bb1.x0 << "," << bb1.x1
          << "] y[" << bb1.y0 << "," << bb1.y1 << "]\n";

    // Explicitly rescale relaxation rate per level (Lagrava et al. 2012, Eq. 24),
    // overriding whatever the generator assigned internally — verifies correctness
    // rather than trusting unconfirmed library behavior.
    //
    // BUG FIX vs the original single-collision-model version: this used to
    // hardcode BGKdynamics on every refined level regardless of
    // collisionModel_runtime -- so an "MRT" run would have been silently
    // BGK on level 1 and above. Now uses the SAME collision model on every
    // level, matching level 0.
    {
        T omega_level = omega_runtime;   // level 0's omega
        for (plint iLevel = 1; iLevel < numLevel; ++iLevel) {
            omega_level = 2.0 * omega_level / (4.0 - omega_level);
            MultiBlockLattice2D<T, DESCRIPTOR>& comp = lattice.getComponent(iLevel);
            defineDynamics(comp, comp.getBoundingBox(), makeBulkDynamics(omega_level, collisionModel_runtime));
            pcout << "Level " << iLevel << " rescaled omega = " << omega_level
                  << " (tau = " << 1.0/omega_level << ", model = " << modelStr << ")\n";
        }
    }

    // airfoilDomain (for initializeMomentumExchange) must be expressed in
    // the FINEST level's own coordinates — scale by 2^interpLevel.
    plint scale = 1 << interpLevel;

    // FIXED: airfoilDomain now derived from the same actual-geometry extent
    // as refineBox above, scaled to the finest level -- the old formula
    // used the fixed param::N_chord (200), which would silently miss
    // force contributions from any part of the grown lattice extending
    // beyond that box.
    Box2D airfoilDomain(refineBox.x0 * scale, refineBox.x1 * scale,
                         refineBox.y0 * scale, refineBox.y1 * scale);

    Array<plint,2> forceIds;
    forceIds[0] = lattice.getComponent(interpLevel).internalStatSubscription().subscribeSum();
    forceIds[1] = lattice.getComponent(interpLevel).internalStatSubscription().subscribeSum();

    // Stamp every level; only the finest level tracks force (see stampAirfoilOnLevel comment above).
    for (plint iLevel = 0; iLevel < numLevel; ++iLevel) {
        stampAirfoilOnLevel(lattice, iLevel, forceIds, iLevel == interpLevel);
    }

    initializeMomentumExchange(lattice.getComponent(interpLevel), airfoilDomain);

    setBoundaryConditions(lattice, AoA_rad_runtime);
    pcout << "Boundary conditions applied\n";

    initializeDomain(lattice, AoA_rad_runtime);
    pcout << "Domain initialised at equilibrium\n";

    pcout << "t=0 check: avgEnergy = " << computeAverageEnergy(lattice.getComponent(0)) << "\n";

    runSimulation(lattice, refineBox, interpLevel, forceIds, AoA_rad_runtime, runTag);

    return 0;
}