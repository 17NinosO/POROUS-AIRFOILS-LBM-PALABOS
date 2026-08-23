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

using namespace plb;    //Palabos namespace
using namespace std;

typedef double T;       //All floating point uses this alias
#define DESCRIPTOR descriptors::MRTD2Q9Descriptor  //2D, 9-velocity lattice

//===================================================================================================
//Simulation Parameters
//===================================================================================================

namespace param {

    //Physical Parameters
    const T Re = 1000.0; //Reynolds number
    const T AoA_deg = 0.0; //Angle of attack [degrees] (default; overridable via argv[1])
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
    const plint numLevel = 3;

    // N_chord: COARSE (level-0) chord length, derived so that the finest
    // level lands exactly on N_chord_fine. This is the constant used by
    // every downstream formula (domain size, geometry normalisation,
    // tau/omega), exactly as in the single-level file. (NOTE: this value
    // depends on BOTH N_chord_fine and numLevel -- with the current
    // settings it is NOT 256; check the printed "N_chord(coarse)" line at
    // runtime rather than trusting a hardcoded comment here.)
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
    //constructor call below (which only takes omega). This is true
    //regardless of collisionModel -- these are documented defaults for a
    //future custom-rate implementation, not active even now that MRT is
    //selected below.
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
    const plint forceLogIter = 10;

    // Cp gets its OWN, LOOSER tolerance -- Cl/Cd are integrated quantities
    // (summed momentum exchange over the whole surface), which smooths out
    // local noise for free. Cp at a single station is a single lattice
    // cell's density sample -- inherently noisier.
    const T cpConvTol = 1e-4;

    // Surface Cp and the convergence check both use genuine ROLLING
    // windows (a circular buffer of the most recent samples), so whenever
    // the run actually stops (early convergence OR hitting maxIter), a
    // valid trailing window is always ready. cpAvgWindow doubles as the
    // convergence-check window.
    const plint cpAvgWindow = 5000;      // rolling window length, iterations
    const plint cpSampleIter = 100;      // Cp sampling cadence within the window

    // Convergence check: every convCheckIter iterations (once past the
    // convMinIter warmup), evaluate relative fluctuation of Cl/Cd AND the
    // max relative fluctuation of surface Cp over the last cpAvgWindow
    // iterations. BOTH must be under tolerance for convRequiredChecks
    // CONSECUTIVE evaluations before declaring convergence and stopping.
    const plint convCheckIter = 1000;
    const plint convMinIter = 20000;
    const int convRequiredChecks = 3;

    // Collision model
    enum CollisionModel { BGK, MRT };
    const CollisionModel collisionModel = MRT;

    //Porosity Parameters
    const bool porous = true;
    const plint overlapWidth = 1;
    const plint behaviorLevel = 0;
    const int n_channels = 4;
    const T channel_half_width = 0.004;   // half-thickness of each channel (chord units)
    const T channel_spacing = 0.02;      // vertical spacing between channel centers (chord units)

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

bool isInsideChannel(plint ix, plint iy) {
    using namespace param;
    if (!porous) return false;
    T y_norm = static_cast<T>(iy - y_foil) / static_cast<T>(N_chord);
    for (int c = 0; c < n_channels; ++c) {
        T offset = (static_cast<T>(c) - (n_channels - 1) / 2.0) * channel_spacing;
        if (std::fabs(y_norm - offset) < channel_half_width) {
            return true;
        }
    }
    return false;
}

bool isInsideAirfoilReal(T x_coarse, T y_coarse) {
    using namespace param;
    T x_norm = (x_coarse - x_foil) / static_cast<T>(N_chord);
    T y_norm = (y_coarse - y_foil) / static_cast<T>(N_chord);
    if (x_norm < 0.0 || x_norm > 1.0) return false;
    T yt = naca0012Thickness(x_norm);
    bool insideBody = (std::fabs(y_norm) < yt);
    return insideBody && !isInsideChannel(static_cast<plint>(std::round(x_coarse)), static_cast<plint>(std::round(y_coarse)));
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

Dynamics<T, DESCRIPTOR>* makeBulkDynamics(T omega_) {
    using namespace param;
    if (collisionModel == BGK)
        return new BGKdynamics<T, DESCRIPTOR>(omega_);
    else
        return new MRTdynamics<T, DESCRIPTOR>(omega_);
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
    T area  = 1.0;

    T lift = -Fx_phys * std::sin(AoA_rad) + Fy_phys * std::cos(AoA_rad);
    T drag =  Fx_phys * std::cos(AoA_rad) + Fy_phys * std::sin(AoA_rad);

    Cl = lift / (q_inf * area);
    Cd = drag / (q_inf * area);
}

//===================================================================================================
//Surface Cp -- WINDOWED MEAN via a genuine ROLLING buffer (circular,
//most-recent-N-samples), not a fixed-endpoint snapshot. Whichever way the
//run ends (early convergence or hitting maxIter), a valid trailing window
//is always ready. Samples at fine-level cells just OUTSIDE the solid
//surface (a few cells offset along y -- AoA is applied via inflow
//direction, not geometry rotation, so the airfoil stays axis-aligned).
//Converts density to Cp via p = rho * cs^2. Uses computeAverageDensity()
//on 1-cell boxes -- MPI-safe regardless of domain decomposition.
//===================================================================================================

class SurfaceCpAccumulator {
public:
    SurfaceCpAccumulator(int nStations_, size_t maxSamples_) : nStations(nStations_),
        maxSamples(maxSamples_), bufUpper(nStations_), bufLower(nStations_), consecutivePasses(0) { }

    void accumulate(MultiGridLattice2D<T, DESCRIPTOR>& lattice, plint interpLevel) {
        using namespace param;
        const plint scale = 1 << interpLevel;
        const T q_inf_lattice = 0.5 * U_lb * U_lb;
        const plint sampleOffsetCells = 3;

        for (int k = 0; k < nStations; ++k) {
            T x_norm = static_cast<T>(k) / static_cast<T>(nStations - 1);
            T yt = naca0012Thickness(x_norm);

            plint ix_fine = static_cast<plint>(std::round((x_foil + x_norm * N_chord) * scale));
            plint iy_upper_fine = static_cast<plint>(std::round((y_foil + yt * N_chord) * scale)) + sampleOffsetCells;
            plint iy_lower_fine = static_cast<plint>(std::round((y_foil - yt * N_chord) * scale)) - sampleOffsetCells;

            T rho_upper = computeAverageDensity(lattice.getComponent(interpLevel),
                                                 Box2D(ix_fine, ix_fine, iy_upper_fine, iy_upper_fine));
            T rho_lower = computeAverageDensity(lattice.getComponent(interpLevel),
                                                 Box2D(ix_fine, ix_fine, iy_lower_fine, iy_lower_fine));

            T Cp_upper = (rho_upper - 1.0) * cs_sq / q_inf_lattice;
            T Cp_lower = (rho_lower - 1.0) * cs_sq / q_inf_lattice;

            bufUpper[k].push_back(Cp_upper);
            if (bufUpper[k].size() > maxSamples) bufUpper[k].pop_front();
            bufLower[k].push_back(Cp_lower);
            if (bufLower[k].size() > maxSamples) bufLower[k].pop_front();
        }
    }

    // Mirrors ForceConvergenceTracker::checkConverged()'s interface. Uses
    // the WORST-CASE (max) relative fluctuation across every station and
    // both surfaces -- a single still-noisy station shouldn't be able to
    // hide behind many quiet ones.
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
        file << "x_over_c,Cp_upper_mean,Cp_upper_std,Cp_lower_mean,Cp_lower_std,n_samples\n";
        for (int k = 0; k < nStations; ++k) {
            T x_norm = static_cast<T>(k) / static_cast<T>(nStations - 1);
            T meanU, stdU, meanL, stdL;
            statsOf(bufUpper[k], meanU, stdU);
            statsOf(bufLower[k], meanL, stdL);
            file << x_norm << "," << meanU << "," << stdU << "," << meanL << "," << stdL << ","
                 << bufUpper[k].size() << "\n";
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
//ValueTracer (prone to false positives: energy averaged over the whole
//domain can plateau even while the near-airfoil flow is still genuinely
//transient). Tracks Cl and Cd directly in a rolling window and requires
//the relative fluctuation of BOTH to stay under tolerance for several
//CONSECUTIVE checks before declaring convergence.
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

// CHANGED vs single-level file: level-0 write is now a much bigger physical
// domain at coarser spacing (still full Lx x Ly, coarse cells). Level>0
// cropping to refineBox*scale is unchanged.
void writeVTK(MultiGridLattice2D<T, DESCRIPTOR>& lattice, plint iter,
              Box2D const& refineBox, Box2D const& innerBox) {
    for (plint iLevel = 0; iLevel < lattice.getNumLevels(); ++iLevel) {
        MultiBlockLattice2D<T, DESCRIPTOR>& c = lattice.getComponent(iLevel);

        Box2D writeBox = c.getBoundingBox();
        if (iLevel == 1) {
            plint scale = 1 << iLevel;
            writeBox = Box2D(refineBox.x0 * scale, refineBox.x1 * scale,
                              refineBox.y0 * scale, refineBox.y1 * scale);
        } else if (iLevel >= 2) {
            // innerBox is already in level-1 coordinates; level iLevel needs
            // it scaled by an additional 2^(iLevel-1) beyond that.
            plint scale = 1 << (iLevel - 1);
            writeBox = Box2D(innerBox.x0 * scale, innerBox.x1 * scale,
                              innerBox.y0 * scale, innerBox.y1 * scale);
        }

        VtkImageOutput2D<T> vtkOut(
            createFileName("vtk_level" + std::to_string(iLevel), iter, 6),
            1.0 / std::pow(2.0, iLevel));
        vtkOut.writeData<float>(*computeVelocityNorm(c, writeBox), "velocityNorm", 1.0);
        vtkOut.writeData<2, float>(*computeVelocity(c, writeBox), "velocity", 1.0);
        vtkOut.writeData<float>(*computeDensity(c, writeBox), "density", 1.0);
    }
}

//===================================================================================================
//Simulation Loop -- convergence now requires BOTH Cl/Cd AND surface Cp to
//be stable (replaces the old domain-averaged-energy ValueTracer). VTK is
//written EXACTLY ONCE, at whichever iteration the loop actually ends on
//(early convergence or maxIter) -- not periodically.
//===================================================================================================
void runSimulation(MultiGridLattice2D<T, DESCRIPTOR>& lattice, Box2D const& refineBox,
                    Box2D const& innerBox, plint interpLevel, Array<plint,2> forceIds, T AoA_rad)
{
    using namespace param;
    std::ofstream forceFile("forces.txt");
    forceFile << "iter,Cl,Cd\n";

    // Both rolling windows sample continuously from iteration 0 (Cl/Cd) or
    // from convMinIter (Cp, gated to avoid wasting MPI-reduction calls on
    // early transient data that would just get flushed out later anyway).
    SurfaceCpAccumulator cpAccum(100, static_cast<size_t>(cpAvgWindow / cpSampleIter));
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
            // Evaluate BOTH unconditionally every time -- short-circuiting
            // on `forceOK && cpAccum.checkConverged(...)` would skip the Cp
            // check whenever forceOK is false, leaving its consecutive-pass
            // streak stale instead of correctly updating/resetting.
            bool forceOK = convTracker.checkConverged(convTol, convRequiredChecks);
            bool cpOK = cpAccum.checkConverged(cpConvTol, convRequiredChecks);

            if (forceOK && cpOK) {
                pcout << "Converged (forces AND Cp both stable) at iter=" << iT
                      << " -- stopping early.\n";
                finalIter = iT;
                convergedEarly = true;
                writeVTK(lattice, iT, refineBox, innerBox);
                break;
            }
        }

        if (iT == maxIter) {
            writeVTK(lattice, iT, refineBox, innerBox);
        }

        lattice.collideAndStream();
    }

    forceFile.close();
    pcout << "Simulation complete at iter=" << finalIter
          << (convergedEarly ? " (converged early)" : " (hit maxIter cap, did not converge)")
          << ". Forces saved to forces.txt\n";

    cpAccum.writeToFile("surface_cp.csv");
}

//===================================================================================================
//Main
//===================================================================================================

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);

    // Force every std::cout write to flush IMMEDIATELY -- default C++
    // stdio block-buffers output whenever it's not connected to an
    // interactive terminal (piped through tee, redirected to a file, or
    // run under PBS). Without this, the new convergence-check output
    // below can sit invisible in a buffer for a long time even though
    // it's genuinely already been printed.
    std::cout << std::unitbuf;

    T AoA_deg_runtime = (argc > 1) ? std::atof(argv[1]) : param::AoA_deg;
    T AoA_rad_runtime = AoA_deg_runtime * M_PI / 180.0;

    pcout << "============================================\n";
    pcout << "   NACA 0012 LBM Simulation - MULTIGRID    \n";
    pcout << "============================================\n";
    pcout << "Re          = " << param::Re       << "\n";
    pcout << "AoA         = " << AoA_deg_runtime  << " deg\n";
    pcout << "numLevel    = " << param::numLevel  << "\n";
    pcout << "N_chord(coarse) = " << param::N_chord << " cells\n";
    pcout << "N_chord(fine, effective) = " << param::N_chord_fine << " cells\n";
    pcout << "Domain(coarse)  = " << param::Lx << " x " << param::Ly << " cells\n";
    pcout << "tau (coarse)    = " << param::tau      << "\n";
    pcout << "omega (coarse)  = " << param::omega    << "\n";
    pcout << "U_lb        = " << param::U_lb     << "\n";
    pcout << "Max iters   = " << param::maxIter  << "\n";
    pcout << "Collision   = "
      << (param::collisionModel == param::BGK ? "BGK" : "MRT") << "\n";
    pcout << "============================================\n";

    writeAirfoilGeometry("airfoil_geometry.csv");
    pcout << "Airfoil geometry written to airfoil_geometry.csv\n";

    plint numLevel = param::numLevel, overlapWidth = param::overlapWidth, behaviorLevel = param::behaviorLevel;
    MultiGridManagement2D management(param::Lx, param::Ly, numLevel, overlapWidth, behaviorLevel);

    // Fine region: airfoil (1 chord) + 5 chord wake downstream, 2 chords
    // upstream margin, 2 chords above/below. All in COARSE (level-0) units;
    // management.refine() and writeVTK both apply the 2^level scaling
    // internally/explicitly where needed.
    Box2D refineBox(param::x_foil - 3*param::N_chord, param::x_foil + 12*param::N_chord,
                     param::y_foil - 4*param::N_chord, param::y_foil + 4*param::N_chord);
    management.refine(0, refineBox);

    plint innerScale = 1 << 1;   // level 1's scale factor
    Box2D innerBox((param::x_foil - param::N_chord/10)*innerScale, (param::x_foil + param::N_chord*11/10)*innerScale,
               (param::y_foil - param::N_chord/8)*innerScale, (param::y_foil + param::N_chord/8)*innerScale);
    management.refine(1, innerBox);

    plint xTiles = global::mpi().getSize(), yTiles = 1, interpLevel = numLevel - 1;
    ParallellizeBySquares2D* parallelizer = new ParallellizeBySquares2D(
        management.getBulks(), management.getBoundingBox(interpLevel), xTiles, yTiles);
    management.parallelize(parallelizer);

    MultiGridLattice2D<T, DESCRIPTOR> lattice =
        *MultiGridGenerator2D<T, DESCRIPTOR>::createRefinedLatticeCubicInterpolationFiltering(
            management, makeBulkDynamics(param::omega), behaviorLevel);
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
    // BUG FIX: this used to hardcode BGKdynamics on every refined level
    // regardless of param::collisionModel -- so with collisionModel=MRT
    // (as set above), levels 1 and 2 were silently running BGK while only
    // level 0 was genuinely MRT. Now uses makeBulkDynamics(), which
    // respects collisionModel consistently on every level.
    {
        std::string modelStr = (param::collisionModel == param::BGK) ? "BGK" : "MRT";
        T omega_level = param::omega;   // level 0's omega
        for (plint iLevel = 1; iLevel < numLevel; ++iLevel) {
            omega_level = 2.0 * omega_level / (4.0 - omega_level);
            MultiBlockLattice2D<T, DESCRIPTOR>& comp = lattice.getComponent(iLevel);
            defineDynamics(comp, comp.getBoundingBox(), makeBulkDynamics(omega_level));
            pcout << "Level " << iLevel << " rescaled omega = " << omega_level
                  << " (tau = " << 1.0/omega_level << ", model = " << modelStr << ")\n";
        }
    }

    // airfoilDomain (for initializeMomentumExchange) must be expressed in
    // the FINEST level's own coordinates — scale by 2^interpLevel.
    plint scale = 1 << interpLevel;

    Box2D airfoilDomain((param::x_foil - 2) * scale, (param::x_foil + param::N_chord + 2) * scale,
                         (param::y_foil - param::N_chord/4) * scale, (param::y_foil + param::N_chord/4) * scale);

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

    runSimulation(lattice, refineBox, innerBox, interpLevel, forceIds, AoA_rad_runtime);

    return 0;
}