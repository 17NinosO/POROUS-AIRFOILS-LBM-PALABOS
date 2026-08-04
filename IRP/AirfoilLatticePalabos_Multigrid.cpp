// ============================================================================
// AirfoilLatticePalabos_Multigrid.cpp
//
// Multigrid (2-level) solver for a parametrised lattice of NACA0012
// airfoils, each independently rotated (theta_ij), acting as a "dynamic
// porous media". Adapted from NACA0012Palabos_Multigrid.cpp -- same LBM
// setup, boundary conditions, and multigrid rescaling pattern as the
// validated single-airfoil solver. The only structural changes are WHAT
// gets stamped as solid and HOW BIG the domain is -- both now come from
// AirfoilLatticeGeometry.h instead of a single fixed airfoil.
//
// GEOMETRY / SOLVER SEPARATION:
//   All geometry (grid layout, per-airfoil angle, margin enforcement) is
//   computed once by geom::buildLattice() in AirfoilLatticeGeometry.h,
//   which has zero Palabos dependency. This file only ever asks that
//   module "is this point inside any airfoil?" -- it never computes
//   airfoil shapes itself. Change the lattice by editing makeLatticeParams()
//   below (mirrors generate_lattice_geometry.cpp); nothing else here needs
//   to change.
//
// FREESTREAM vs PER-AIRFOIL ANGLE -- IMPORTANT DISTINCTION:
//   As in the single-airfoil file, the FREESTREAM direction (AoA_inflow)
//   is applied by rotating the INLET VELOCITY -- the domain stays axis-
//   aligned. Each individual airfoil's angle (theta_ij) is a SEPARATE,
//   independent rotation applied to that airfoil's own geometry inside
//   AirfoilLatticeGeometry.h. Sweep duct AoA and the per-airfoil "porous"
//   angle field independently.
//
// FORCE TRACKING -- v1 SCOPE:
//   This measures ONE combined force over the WHOLE lattice (bulk drag/
//   lift on the array) -- the porous-media analogue of a single airfoil's
//   force. Splitting into per-airfoil forces (one MomentumExchangeBounceBack
//   + one forceIds pair per unit) is a natural follow-on, not implemented
//   here yet -- see the comment above stampLatticeOnLevel().
//
// BUILD (same toolchain/modules as your existing solver on Crescent2):
//   Place AirfoilLatticeGeometry.h in the same folder (or add its folder
//   to the include path) and compile exactly like NACA0012Palabos_Multigrid.cpp,
//   just with this file as the source.
// ============================================================================

#include "palabos2D.h"  //Declarations: class interfaces, types
#include "palabos2D.hh" //Definitions: template implementations

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "complexDynamics/mrtDynamics.h"
#include "boundaryCondition/bounceBackModels.h"
#include "complexDynamics/mrtDynamics.hh"

#include "AirfoilLatticeGeometry.h"

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::MRTD2Q9Descriptor

//===================================================================================================
//Simulation Parameters (unchanged in spirit from NACA0012Palabos_Multigrid.cpp)
//===================================================================================================

namespace param {

    //Physical Parameters
    const T Re = 1000.0; //Reynolds number
    const T AoA_inflow_deg = 0.0; //Freestream/duct angle [deg] (default; overridable via argv[1])
    const T AoA_inflow_rad = AoA_inflow_deg * M_PI / 180.0;

    //---------------------------------------------------------------------
    // MULTIGRID RESOLUTION SETUP (identical convention to the single-airfoil file)
    //---------------------------------------------------------------------
    const plint N_chord_fine = 400;
    const plint numLevel = 2;
    const plint N_chord = N_chord_fine / (1 << (numLevel - 1));

    //Lattice velocity
    const T U_lb = 0.05;

    //Derived LBM Parameters (computed at the COARSE/behaviorLevel resolution)
    const T cs_sq = 1.0 / 3.0;
    const T nu_lb = U_lb * N_chord / Re;
    const T tau = 0.5 + nu_lb / cs_sq;
    const T omega = 1.0 / tau;

    //Simulation Control Parameters
    const plint maxIter = 100000;
    const plint outIter = 1000;
    const T convTol = 1e-6;
    const plint forceLogIter = 10;

    enum CollisionModel { BGK, MRT };
    const CollisionModel collisionModel = BGK;

    const plint overlapWidth = 1;
    const plint behaviorLevel = 0;

} //namespace param

//===================================================================================================
//Lattice geometry parameters -- feeds straight into geom::buildLattice().
//Keep these in sync with generate_lattice_geometry.cpp if you want the
//standalone plot to match what the solver actually runs.
//===================================================================================================

geom::LatticeParams makeLatticeParams() {
    geom::LatticeParams p;
    p.nx = 6;
    p.ny = 4;
    p.N_chord = param::N_chord; // COARSE chord -- shared with the solver, not N_chord_fine

    // X: 1/4-chord GAP between one airfoil's TE and the next's LE.
    //    pitch = chord (1.0) + gap (0.25)
    p.pitch_x_chords = 1.25;
    // Y: 1/2-chord GAP between rows, measured against NACA0012's max
    //    thickness (0.12c) since that's the airfoil's own extent in y.
    //    pitch = max_thickness (0.12) + gap (0.5)
    p.pitch_y_chords = 0.62;

    p.pivot_chords = 0.25; // quarter-chord rotation pivot

    // Margins now scale with the WHOLE ARRAY's own footprint, not one
    // chord -- treats the lattice as a single bluff body. This OVERWRITES
    // margin_inlet_chords / margin_outlet_chords / margin_wall_chords;
    // tune the multipliers (array widths for inlet/outlet, array heights
    // for wall) rather than the margins directly. Must be called before
    // geom::buildLattice().
    // inlet_widths=3.0 matches the original single-airfoil code's
    // "3 chords upstream" convention, now applied to the array's own
    // width instead of one chord.
    // --- boundary clearance, scaled to the WHOLE ARRAY's own footprint
    // (treating the lattice as one bluff body) rather than a single
    // chord. Keeps the original solver's asymmetry (inlet closer, outlet
    // further for wake development) but at multipliers that stay
    // compute-feasible when the characteristic length is the whole array
    // (~7.9c wide) rather than one chord -- the literal 8/27/8 ratio
    // would give a 670M-cell domain at N_chord=200. ---
    geom::applyBodyScaledMargins(p, /*inlet_widths=*/2.0, /*outlet_widths=*/4.0, /*wall_heights=*/2.0);

    // Placeholder angle field -- replace with a lookup table or an
    // optimiser's output when ready. Keep identical to
    // generate_lattice_geometry.cpp if you want the two to match.
    p.theta_deg_fn = [](int i, int j) {
        return 15.0 * std::sin(0.6 * j) * std::cos(0.4 * i);
    };

    return p;
}

//===================================================================================================
//Flags (unchanged)
//===================================================================================================

namespace cellFlag {
    const int FLUID = 0;
    const int SOLID = 1;
    const int INLET = 2;
    const int OUTLET = 3;
    const int WALL = 4;
}

//===================================================================================================
//Stamping -- same role as SetAirfoilSolidOnLevel in the single-airfoil
//file, but tests against the WHOLE lattice via geom::isInsideLattice()
//instead of one fixed airfoil.
//===================================================================================================

class SetLatticeSolidOnLevel : public BoxProcessingFunctional2D_S<int> {
public:
    SetLatticeSolidOnLevel(plint scale_, const geom::LatticeGeometry& g_,
                            plint N_chord_, T pivot_chords_)
        : scale(scale_), g(g_), N_chord(N_chord_), pivot_chords(pivot_chords_) { }

    virtual void process(Box2D domain, ScalarField2D<int>& flags) {
        Dot2D offset = flags.getLocation();
        for (plint ix = domain.x0; ix <= domain.x1; ++ix) {
            for (plint iy = domain.y0; iy <= domain.y1; ++iy) {
                plint gx = ix + offset.x;
                plint gy = iy + offset.y;
                T x_coarse = static_cast<T>(gx) / scale;
                T y_coarse = static_cast<T>(gy) / scale;
                if (geom::isInsideLattice(x_coarse, y_coarse, g.units, N_chord, pivot_chords)) {
                    flags.get(ix, iy) = cellFlag::SOLID;
                }
            }
        }
    }

    virtual SetLatticeSolidOnLevel* clone() const {
        return new SetLatticeSolidOnLevel(*this);
    }

    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
    }

private:
    plint scale;
    geom::LatticeGeometry g;
    plint N_chord;
    T pivot_chords;
};

// Same trackForce convention as the single-airfoil file: only the finest
// level gets MomentumExchangeBounceBack (force-accumulating); every other
// level gets plain BounceBack, so the array's footprint stays geometrically
// consistent under the fine patch without double-counting force into
// forceIds from more than one level.
void stampLatticeOnLevel(MultiGridLattice2D<T, DESCRIPTOR>& lattice, plint iLevel,
                          const geom::LatticeGeometry& g, T pivot_chords,
                          Array<plint,2> forceIds, bool trackForce) {
    MultiBlockLattice2D<T, DESCRIPTOR>& comp = lattice.getComponent(iLevel);
    Box2D bb = comp.getBoundingBox();
    plint scale = 1 << iLevel;

    MultiScalarField2D<int> flags(comp);
    applyProcessingFunctional(
        new SetLatticeSolidOnLevel(scale, g, param::N_chord, pivot_chords), bb, flags);

    if (trackForce) {
        defineDynamics(comp, flags, new MomentumExchangeBounceBack<T, DESCRIPTOR>(forceIds), cellFlag::SOLID);
    } else {
        defineDynamics(comp, flags, new BounceBack<T, DESCRIPTOR>(), cellFlag::SOLID);
    }
    pcout << "Lattice stamped on level " << iLevel
          << (trackForce ? " (force-tracked)" : " (plain, no force tracking)") << "\n";
}

//===================================================================================================
//Dynamics factory (unchanged)
//===================================================================================================

Dynamics<T, DESCRIPTOR>* makeBulkDynamics() {
    using namespace param;
    if (collisionModel == BGK)
        return new BGKdynamics<T, DESCRIPTOR>(omega);
    else
        return new MRTdynamics<T, DESCRIPTOR>(omega);
}

//===================================================================================================
//Initialisation (unchanged; driven by AoA_inflow_rad, NOT any theta_ij)
//===================================================================================================

void initializeDomain(MultiGridLattice2D<T, DESCRIPTOR>& lattice, T AoA_inflow_rad) {
    using namespace param;
    const T rho_0 = 1.0;
    Array<T,2> u_inf(U_lb * std::cos(AoA_inflow_rad), U_lb * std::sin(AoA_inflow_rad));

    for (plint iLevel = 0; iLevel < lattice.getNumLevels(); ++iLevel) {
        initializeAtEquilibrium(
            lattice.getComponent(iLevel),
            lattice.getComponent(iLevel).getBoundingBox(),
            rho_0, u_inf);
    }
    lattice.initialize();
}

//===================================================================================================
//Boundary Conditions -- identical pattern to the single-airfoil file, now
//sized off the lattice's own DERIVED Lx, Ly instead of a fixed formula.
//===================================================================================================

void setBoundaryConditions(MultiGridLattice2D<T, DESCRIPTOR>& lattice, plint Lx, plint Ly, T AoA_inflow_rad) {
    using namespace param;
    Array<T, 2> u_inf(U_lb * std::cos(AoA_inflow_rad), U_lb * std::sin(AoA_inflow_rad));

    MultiBlockLattice2D<T, DESCRIPTOR>& coarse = lattice.getComponent(0);

    OnLatticeBoundaryCondition2D<T, DESCRIPTOR>* bc =
        createLocalBoundaryCondition2D<T, DESCRIPTOR>();

    // Inlet (left) -- velocity
    bc->addVelocityBoundary0N(Box2D(0, 0, 1, Ly-2), coarse);
    setBoundaryVelocity(coarse, Box2D(0, 0, 1, Ly-2), u_inf);

    // Outlet (right) -- pressure
    bc->addPressureBoundary0P(Box2D(Lx-1, Lx-1, 1, Ly-2), coarse);
    setBoundaryDensity(coarse, Box2D(Lx-1, Lx-1, 1, Ly-2), (T)1.0);

    // Top and bottom -- velocity, matching free-stream
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
//Force Computation -- decomposed against the FREESTREAM direction
//(AoA_inflow_rad), independent of any individual airfoil's theta_ij.
//Reports the bulk force on the WHOLE array; "per-airfoil average" divides
//by n_foils purely as a sanity-check quantity against your single-airfoil
//baseline results, NOT a physically rigorous per-unit load.
//===================================================================================================

void computeCoefficients(T Fx, T Fy, T& Cl_total, T& Cd_total, T& Cl_avg, T& Cd_avg,
                          T AoA_inflow_rad, int n_foils) {
    using namespace param;
    T dx = 1.0 / (T)N_chord_fine;
    T dt = U_lb / (T)N_chord_fine;
    T scale = dx*dx*dx / (dt*dt);

    T Fx_phys = 2.0 * Fx * scale;
    T Fy_phys = 2.0 * Fy * scale;

    T q_inf = 0.5 * 1.0 * 1.0;
    T area  = 1.0; // one chord, per unit span

    T lift = -Fx_phys * std::sin(AoA_inflow_rad) + Fy_phys * std::cos(AoA_inflow_rad);
    T drag =  Fx_phys * std::cos(AoA_inflow_rad) + Fy_phys * std::sin(AoA_inflow_rad);

    Cl_total = lift / (q_inf * area);
    Cd_total = drag / (q_inf * area);
    Cl_avg = Cl_total / static_cast<T>(n_foils);
    Cd_avg = Cd_total / static_cast<T>(n_foils);
}

//===================================================================================================
//VTK output (unchanged pattern)
//===================================================================================================

void writeVTK(MultiGridLattice2D<T, DESCRIPTOR>& lattice, plint iter, Box2D const& refineBox) {
    for (plint iLevel = 0; iLevel < lattice.getNumLevels(); ++iLevel) {
        MultiBlockLattice2D<T, DESCRIPTOR>& c = lattice.getComponent(iLevel);

        Box2D writeBox = c.getBoundingBox();
        if (iLevel > 0) {
            plint scale = 1 << iLevel;
            writeBox = Box2D(refineBox.x0 * scale, refineBox.x1 * scale,
                              refineBox.y0 * scale, refineBox.y1 * scale);
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
//Simulation Loop (unchanged pattern)
//===================================================================================================

void runSimulation(MultiGridLattice2D<T, DESCRIPTOR>& lattice, Box2D const& refineBox,
                    plint interpLevel, Array<plint,2> forceIds, T AoA_inflow_rad, int n_foils)
{
    using namespace param;
    std::ofstream forceFile("forces.txt");
    forceFile << "iter,Cl_total,Cd_total,Cl_avg_per_foil,Cd_avg_per_foil\n";

    util::ValueTracer<T> convergence(param::U_lb, (T)param::N_chord, param::convTol);

    pcout << "Starting simulation: " << maxIter << " max iterations\n";

    for (plint iT = 0; iT <= maxIter; ++iT) {

        if (iT % forceLogIter == 0) {
            T Fx = lattice.getComponent(interpLevel).getInternalStatistics().getSum(forceIds[0]);
            T Fy = lattice.getComponent(interpLevel).getInternalStatistics().getSum(forceIds[1]);

            const T subStepsPerCoarseIter = static_cast<T>(1 << interpLevel);
            Fx /= subStepsPerCoarseIter;
            Fy /= subStepsPerCoarseIter;

            T Cl_total, Cd_total, Cl_avg, Cd_avg;
            computeCoefficients(Fx, Fy, Cl_total, Cd_total, Cl_avg, Cd_avg, AoA_inflow_rad, n_foils);
            forceFile << iT << "," << Cl_total << "," << Cd_total << ","
                       << Cl_avg << "," << Cd_avg << "\n";
            forceFile.flush();
        }

        if (iT % outIter == 0) {
            T avgEnergy = computeAverageEnergy(lattice.getComponent(0));
            pcout << "iter=" << iT << " E=" << avgEnergy << "\n";
            convergence.takeValue(avgEnergy, true);
            if (convergence.hasConverged()) {
                pcout << "Converged at iter=" << iT << "\n";
                writeVTK(lattice, iT, refineBox);
                break;
            }
        }

        if (iT % 10000 == 0 || iT == maxIter) {
            writeVTK(lattice, iT, refineBox);
        }

        lattice.collideAndStream();
    }

    forceFile.close();
    pcout << "Simulation complete. Forces saved to forces.txt\n";
}

//===================================================================================================
//Main
//===================================================================================================

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);

    T AoA_inflow_deg_runtime = (argc > 1) ? std::atof(argv[1]) : param::AoA_inflow_deg;
    T AoA_inflow_rad_runtime = AoA_inflow_deg_runtime * M_PI / 180.0;

    geom::LatticeParams latticeParams = makeLatticeParams();
    geom::LatticeGeometry g = geom::buildLattice(latticeParams);
    int n_foils = latticeParams.nx * latticeParams.ny;

    plint Lx = g.Lx;
    plint Ly = g.Ly;

    // Overlap check -- tighter pitch trades away the conservative
    // bounding-circle safety margin, so verify against actual rotated
    // surfaces before running anything.
    {
        std::vector<std::pair<int,int>> overlaps =
            geom::findOverlaps(g, param::N_chord, latticeParams.pivot_chords);
        if (!overlaps.empty()) {
            pcout << "WARNING: " << overlaps.size()
                  << " airfoil pair(s) geometrically overlap at this pitch/angle "
                  << "combination -- results will not be physical. Increase pitch "
                  << "or reduce the angle field.\n";
            for (auto& pr : overlaps) {
                pcout << "  units[" << pr.first << "] x units[" << pr.second << "]\n";
            }
        }
    }

    pcout << "============================================\n";
    pcout << "   NACA 0012 AIRFOIL LATTICE - MULTIGRID    \n";
    pcout << "============================================\n";
    pcout << "Re              = " << param::Re       << "\n";
    pcout << "AoA (inflow)    = " << AoA_inflow_deg_runtime << " deg\n";
    pcout << "Lattice         = " << latticeParams.nx << " x " << latticeParams.ny
          << " = " << n_foils << " airfoils\n";
    pcout << "numLevel        = " << param::numLevel  << "\n";
    pcout << "N_chord(coarse) = " << param::N_chord << " cells\n";
    pcout << "N_chord(fine, effective) = " << param::N_chord_fine << " cells\n";
    pcout << "Domain(coarse)  = " << Lx << " x " << Ly << " cells\n";
    pcout << "tau (coarse)    = " << param::tau      << "\n";
    pcout << "omega (coarse)  = " << param::omega    << "\n";
    pcout << "U_lb            = " << param::U_lb     << "\n";
    pcout << "Max iters       = " << param::maxIter  << "\n";
    pcout << "Collision       = "
      << (param::collisionModel == param::BGK ? "BGK" : "MRT") << "\n";
    pcout << "============================================\n";

    geom::writeLatticeCenters(g, "lattice_centers.csv");
    geom::writeLatticeSurfaces(g, param::N_chord, "lattice_surfaces.csv");
    pcout << "Lattice geometry written to lattice_centers.csv / lattice_surfaces.csv\n";

    plint numLevel = param::numLevel, overlapWidth = param::overlapWidth, behaviorLevel = param::behaviorLevel;
    MultiGridManagement2D management(Lx, Ly, numLevel, overlapWidth, behaviorLevel);

    // Fine region: hugs the array's OWN extent, with a small fixed padding
    // -- deliberately decoupled from margin size. Now that margins scale
    // with the whole array (potentially large), tying the refine box to
    // margin size would waste huge amounts of fine-grid compute on empty
    // buffer space; tying it to the array's own extent instead keeps level
    // 1 tight in both x and y no matter how generous the physical margins
    // are. This also means the refine box can never reach the domain edge
    // (the margins guarantee real coarse-grid space beyond it), so the
    // earlier boundary-touching multigrid segfault can't recur here.
    T C = static_cast<T>(param::N_chord);
    T R = std::max(latticeParams.pivot_chords, 1.0 - latticeParams.pivot_chords) + 0.06;

    // Reach of the array itself (chord units), independent of margins:
    T xIn_chords  = latticeParams.margin_inlet_chords; // leftmost reach = margin_inlet, by construction (see buildLattice)
    T xOut_chords = latticeParams.margin_inlet_chords + 2.0 * R
                    + (latticeParams.nx - 1) * latticeParams.pitch_x_chords;
    T yBot_chords = latticeParams.margin_wall_chords;
    T yTop_chords = latticeParams.margin_wall_chords + 2.0 * R
                    + (latticeParams.ny - 1) * latticeParams.pitch_y_chords;

    const T padNear = 1.5; // chords, tight fixed padding on every side
    const T padWake = 3.0; // chords, small extra allowance for the near wake

    plint xmin = static_cast<plint>((xIn_chords  - padNear) * C);
    plint xmax = static_cast<plint>((xOut_chords + padNear + padWake) * C);
    plint ymin = static_cast<plint>((yBot_chords - padNear) * C);
    plint ymax = static_cast<plint>((yTop_chords + padNear) * C);

    Box2D refineBox(std::max(xmin, (plint)1), std::min(xmax, Lx - 2),
                     std::max(ymin, (plint)1), std::min(ymax, Ly - 2));
    management.refine(0, refineBox);

    plint xTiles = global::mpi().getSize(), yTiles = 1, interpLevel = numLevel - 1;
    ParallellizeBySquares2D* parallelizer = new ParallellizeBySquares2D(
        management.getBulks(), management.getBoundingBox(interpLevel), xTiles, yTiles);
    management.parallelize(parallelizer);

    MultiGridLattice2D<T, DESCRIPTOR> lattice =
        *MultiGridGenerator2D<T, DESCRIPTOR>::createRefinedLatticeCubicInterpolationFiltering(
            management, makeBulkDynamics(), behaviorLevel);
    pcout << "Lattice built (numLevel=" << numLevel << ")\n";

    Box2D bb0 = lattice.getComponent(0).getBoundingBox();
    pcout << "Level 0 bbox: x[" << bb0.x0 << "," << bb0.x1 << "] y[" << bb0.y0 << "," << bb0.y1 << "]\n";

    Box2D bb1 = lattice.getComponent(interpLevel).getBoundingBox();
    pcout << "Level " << interpLevel << " bbox: x[" << bb1.x0 << "," << bb1.x1
          << "] y[" << bb1.y0 << "," << bb1.y1 << "]\n";

    // Explicitly rescale relaxation rate per level (Lagrava et al. 2012, Eq. 24),
    // same as the single-airfoil file.
    {
        T omega_level = param::omega;
        for (plint iLevel = 1; iLevel < numLevel; ++iLevel) {
            omega_level = 2.0 * omega_level / (4.0 - omega_level);
            MultiBlockLattice2D<T, DESCRIPTOR>& comp = lattice.getComponent(iLevel);
            defineDynamics(comp, comp.getBoundingBox(), new BGKdynamics<T, DESCRIPTOR>(omega_level));
            pcout << "Level " << iLevel << " rescaled omega = " << omega_level
                  << " (tau = " << 1.0/omega_level << ")\n";
        }
    }

    Array<plint,2> forceIds;
    forceIds[0] = lattice.getComponent(interpLevel).internalStatSubscription().subscribeSum();
    forceIds[1] = lattice.getComponent(interpLevel).internalStatSubscription().subscribeSum();

    // Stamp every level; only the finest level tracks force (see
    // stampLatticeOnLevel comment above).
    for (plint iLevel = 0; iLevel < numLevel; ++iLevel) {
        stampLatticeOnLevel(lattice, iLevel, g, latticeParams.pivot_chords, forceIds, iLevel == interpLevel);
    }

    // Momentum-exchange domain: the whole refine box, expressed in the
    // finest level's own coordinates -- same role as airfoilDomain in the
    // single-airfoil file, just sized to cover the full lattice footprint.
    plint scale = 1 << interpLevel;
    Box2D latticeDomain(refineBox.x0 * scale, refineBox.x1 * scale,
                         refineBox.y0 * scale, refineBox.y1 * scale);
    initializeMomentumExchange(lattice.getComponent(interpLevel), latticeDomain);

    setBoundaryConditions(lattice, Lx, Ly, AoA_inflow_rad_runtime);
    pcout << "Boundary conditions applied\n";

    initializeDomain(lattice, AoA_inflow_rad_runtime);
    pcout << "Domain initialised at equilibrium\n";

    pcout << "t=0 check: avgEnergy = " << computeAverageEnergy(lattice.getComponent(0)) << "\n";

    runSimulation(lattice, refineBox, interpLevel, forceIds, AoA_inflow_rad_runtime, n_foils);

    return 0;
}