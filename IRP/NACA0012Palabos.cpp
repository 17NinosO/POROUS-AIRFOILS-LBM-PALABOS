#include "palabos2D.h"  //Declarations:class interfaces, types
#include "palabos2D.hh" //Definitions: template implementations

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>

#include "complexDynamics/mrtDynamics.h"
#include "boundaryCondition/bounceBackModels.h"
#include "complexDynamics/mrtDynamics.hh"

using namespace plb;    //Palabos namespace
using namespace std;

// Palabos is a template library (e.g. MultiGridLattice2D<T, DESCRIPTOR>).
//C++ requires template definitions to be visible at compile time, not just declarations.
//Palabos ships .hh files containing full implementations, separate from the .h declaration headers.

//Precisiion type and lattice descriptor.

typedef double T;       //All floating point uses this alias
#define DESCRIPTOR descriptors::MRTD2Q9Descriptor  //2D, 9-velocity lattice

//D2Q9 means that there are 2 dimensions and 9 velocity directions.
//LBM has a finite set of discrete velocity vectors.
//D2Q9 has one rest velocity, 4 axis-aligned, and 4 diagonal velocities.
//These are enough to recover the NS equations to a second order.
//This is a standard choice for a 2D incompressible flow around an airfoil.

//The LBM unit system.
//LBM doesn't work in SI unites but they are expressed in lattice units.
//Lattice spacing Delta_x = 1
//Time step Delta_t = 1
//Lattice speed of sound c_s = 1/sqrt(3)

//A conversion is defined between the physical and lattice quantities through the Re number.

//Re = U_lb * N_chord / v_lb

//N_chord represents the chord length in lattice cells.

//===================================================================================================
//Simulation Parameters
//===================================================================================================

namespace param {

    //Physical Parameters
    const T Re = 1000.0; //Reynolds number
    const T AoA_deg = 0.0; //Angle of attack [degrees]
    const T AoA_rad = AoA_deg * M_PI / 180.0; //Angle of attack [radians]

    //Lattice Resolution
    //N_chord: number of lattice cells spanning the chord
    //Higher = more accurate but more computationally expensive
    //100 is a good starting point; but use 50 for a quick test run.
    const plint N_chord = 1024;

    //Lattice velocity
    //Must satisfy Ma__lb = U_lb / c_s << 1 for incompressible flow.
    //c_s = 1/sqrt(3) ~ 0.577 in lattice units.
    const T U_lb = 0.05;

    //Derived LBM Parameters
    const T cs_sq = 1.0 / 3.0;
    const T nu_lb = U_lb * N_chord / Re; //Kinematic viscosity in lattice units.
    const T tau = 0.5 + nu_lb / cs_sq; //Relaxation time for BGK collision operator.
    const T omega = 1.0 / tau; //Relaxation rate for BGK collision operator.

    //MRT relaxation rates (D2Q9)
    //Omega controls viscosity - fixed by Re
    //Other rates are free parameters tuned for stability
    //conserved moments = 0, non-hydro moments 0<s<2
    const T s_1 = 1.4;
    const T s_2 = 1.4;
    const T s_4 = 1.2;
    const T s_6 = 1.2;
    //s_7 = s_8 = omega (shear stress - sets viscosity)
    

    //Domain Size [Lattice Unites]
    //Simple Rectangular domain: inlet left, outlet right, slip walls top and bottom.
    //8 chord lengths upstream and 16 chord lengths downstream of the airfoil.
    const plint Lx = 36 * N_chord; // Domain Width
    const plint Ly = 16 * N_chord; // Domain Height 
    const plint x_foil = 8 * N_chord; // LE x position from the inlet.
    const plint y_foil = Ly / 2; // LE y centred vertically.

    //Simulation Control Parameters
    const plint maxIter = 200000; //Maximum number of iterations.
    const plint outIter = 1000; //Output and Logging Interval
    const T convTol = 1e-6; //Velocity convergence tolerance.
    const plint forceLogIter = 10; //Force/Cl/Cd logging interval

    // Collision model
    enum CollisionModel { BGK, MRT };
    const CollisionModel collisionModel = BGK;   //change BGK/MRT

    //Porosity Parameters
    const bool porous = false;
    // Grid refinement (Stage 1: coarse + one fine level, level-0 geometry only)
    const plint numLevel     = 1;
    const plint overlapWidth = 1;
    const plint behaviorLevel = 0;
    const T pore_width = 0.02;
    const int n_pores = 4;
    const T pore_centre_x = 0.30;
    const T pore_spacing = 0.08;

} //namespace param

//===================================================================================================
//Airfoil Geometry
//===================================================================================================

//This section generates the airfoil surface coordinates analytically.
//And provides an inside/outside test so it can flag solid lattice cells later.

//Angle of attack will be applied by changing the inlet velocity direction, not by rotating the airfoil geometry.

// NACA GEOMETRY

//Half thickness of NACA 0012 at normalised chord position x [0,1]
//Returns y_t normalised by chord (multiply by N_chord to get lattice units).
T naca0012Thickness(T x) {
    if (x<0.0) x=0.0;
    if (x>1.0) x=1.0;

    //NACA 4 digit coefficients
    // a4 = 0.1015 gives an open trailing edge (y_t ~ 0.00126 at x=1)
    //difference is < 1 lattice cell
    return 5.0 * 0.12 * (
        0.2969 * std::sqrt(x)
        - 0.1260 * x
        -0.3516 * x*x
        + 0.2843 * x*x*x
        - 0.1015 * x*x*x*x
    );
}

//Pore Geometry Test

bool isInsidePore(plint ix, plint iy) {
    using namespace param;
    
    //If porosity is switched off, no cell is ever inside a pore
    if (!porous) return false;

    //Map the cell to normalised chord coordinates
    T x_norm = static_cast<T>(ix - x_foil) / static_cast<T>(N_chord);

    //Loop over each pore in the cluster
    for (int p = 0; p< n_pores; ++p) {

        T offset    =(static_cast<T>(p) - (n_pores - 1) / 2.0) * pore_spacing;
        T pore_x    = pore_centre_x + offset;

        if (std::fabs(x_norm - pore_x) < pore_width / 2.0){
            return true;
        }
    }

    return false;
}

//Inside/Outside test for the airfoil geometry.
//For every lattice cell (iX, iY) is inside the NACA solid
//The airfoil leading edge sits at (param::x_foil, param::y_foil).

bool isInsideAirfoilReal(T x_coarse, T y_coarse) {
    using namespace param;
    T x_norm = (x_coarse - x_foil) / static_cast<T>(N_chord);
    T y_norm = (y_coarse - y_foil) / static_cast<T>(N_chord);
    if (x_norm < 0.0 || x_norm > 1.0) return false;
    T yt = naca0012Thickness(x_norm);
    bool insideBody = (std::fabs(y_norm) < yt);
    return insideBody && !isInsidePore(static_cast<plint>(std::round(x_coarse)), static_cast<plint>(std::round(y_coarse)));
}

namespace cellFlag {
    const int FLUID = 0;    //Regular fluid cell
    const int SOLID = 1;    //Airfoil interior
    const int INLET = 2;    //Left boundary
    const int OUTLET = 3;   //Right boundary
    const int WALL = 4;     //Top and bottom boundaries
}

//Stamps SOLID flags for cells inside the airfoil, for one grid level.
//Runs as a per-block functional (not a raw MultiScalarField2D::get() loop)
//since get() on a field sharing a multi-block's real block structure
//pays a per-call block-lookup cost that is ruinous across millions of cells.
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

void stampAirfoilOnLevel(MultiGridLattice2D<T, DESCRIPTOR>& lattice, plint iLevel, Array<plint,2> forceIds) {
    MultiBlockLattice2D<T, DESCRIPTOR>& comp = lattice.getComponent(iLevel);
    Box2D bb = comp.getBoundingBox();
    plint scale = 1 << iLevel;

    MultiScalarField2D<int> flags(comp);
    applyProcessingFunctional(new SetAirfoilSolidOnLevel(scale), bb, flags);
    defineDynamics(comp, flags, new MomentumExchangeBounceBack<T, DESCRIPTOR>(forceIds), cellFlag::SOLID);
    pcout << "Airfoil stamped on level " << iLevel << "\n";
}

//===================================================================================================
//Surface Point Generator (For output and post-processing)
//===================================================================================================

//This section writes the NACA coordinates to a CSV file
//Use cosinne spacing to cluster points near the leading and trailing edges where curvature is highest.

void writeAirfoilGeometry(const std::string& filename, int nPoints=200) {
    using namespace param;

    std::ofstream file(filename);
    file << "x_upper, y_upper, x_lower, y_lower\n";

    for (int i = 0; i < nPoints; ++i) {
        //cosine spacing:theta sweeps 0 to pi, giving x_norm 0-1
        T theta = M_PI * static_cast<T>(i) / static_cast<T>(nPoints - 1);
        T x_norm = 0.5 * (1.0 - std::cos(theta));

        T yt = naca0012Thickness(x_norm);

        //Scale to lattice units

        T x_lat = x_foil + x_norm * N_chord;
        T y_upper = y_foil + yt * N_chord;
        T y_lower = y_foil - yt * N_chord;

        file << x_lat << "," << y_upper << ","
            << x_lat << "," << y_lower << "\n";

    }
}

//===================================================================================================
//Geometry Flags and Lattice Construction
//===================================================================================================

//This section builds the flag matrix which is a scalar field that labels every cell as a fluid, solid, inlet, outlet or wall.
//Creating the MultiGridLattice and stamping the correct collision dynamics onto each cell type.

// Returns a fresh dynamics object of the chosen type.
// Called wherever the code needs "the bulk fluid dynamics".
Dynamics<T, DESCRIPTOR>* makeBulkDynamics() {
    using namespace param;
    if (collisionModel == BGK)
        return new BGKdynamics<T, DESCRIPTOR>(omega);
    else
        return new MRTdynamics<T, DESCRIPTOR>(omega);
}

//Creating the lattice and assigning dynamics
//Create the MultiGridLattice2D, and assign dynamics from the flag matrix.
//Returns a raw pointer

//===================================================================================================
//Initialisation
//===================================================================================================

void initializeDomain(MultiGridLattice2D<T, DESCRIPTOR>& lattice) {
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
//Boundary Conditions
//===================================================================================================
void setBoundaryConditions(MultiGridLattice2D<T, DESCRIPTOR>& lattice) {
    using namespace param;
    Array<T, 2> u_inf(U_lb * std::cos(AoA_rad), U_lb * std::sin(AoA_rad));

    // Domain-edge boundaries live on the coarse level only.
    MultiBlockLattice2D<T, DESCRIPTOR>& coarse = lattice.getComponent(0);

    OnLatticeBoundaryCondition2D<T, DESCRIPTOR>* bc =
        createLocalBoundaryCondition2D<T, DESCRIPTOR>();

    bc->addVelocityBoundary0N(Box2D(0, 0, 1, Ly-2), coarse);
    setBoundaryVelocity(coarse, Box2D(0, 0, 1, Ly-2), u_inf);
    bc->addPressureBoundary0P(Box2D(Lx-1, Lx-1, 1, Ly-2), coarse);
    setBoundaryDensity(coarse, Box2D(Lx-1, Lx-1, 1, Ly-2), (T)1.0);
    bc->addVelocityBoundary1P(Box2D(1, Lx-2, Ly-1, Ly-1), coarse);
    setBoundaryVelocity(coarse, Box2D(1, Lx-2, Ly-1, Ly-1), u_inf);
    bc->addVelocityBoundary1N(Box2D(1, Lx-2, 0, 0), coarse);
    setBoundaryVelocity(coarse, Box2D(1, Lx-2, 0, 0), u_inf);

    bc->addPressureBoundary1P(Box2D(1, Lx-2, Ly-1, Ly-1), coarse);
    setBoundaryDensity(coarse, Box2D(1, Lx-2, Ly-1, Ly-1), (T)1.0);
    bc->addPressureBoundary1N(Box2D(1, Lx-2, 0, 0), coarse);
    setBoundaryDensity(coarse, Box2D(1, Lx-2, 0, 0), (T)1.0);

    defineDynamics(coarse, Box2D(0,0,0,0),       new BounceBack<T,DESCRIPTOR>());
    defineDynamics(coarse, Box2D(0,0,Ly-1,Ly-1), new BounceBack<T,DESCRIPTOR>());
    defineDynamics(coarse, Box2D(Lx-1,Lx-1,0,0), new BounceBack<T,DESCRIPTOR>());
    defineDynamics(coarse, Box2D(Lx-1,Lx-1,Ly-1,Ly-1), new BounceBack<T,DESCRIPTOR>());

    delete bc;
}

//===================================================================================================
//Force Computation
//===================================================================================================

//Computes raw lattice-unit forces Fx, Fy on the airfoil
//using the momentum exchange method
//iterates only over the airfoil bounding box

void computeCoefficients(T Fx, T Fy, T& Cl, T& Cd) {
    using namespace param;

    // Lattice-unit conversion (matches Gabriel's IncomprFlowParam convention)
    T dx = 1.0 / (T)N_chord;          // grid spacing, chord normalised to 1
    T dt = U_lb / (T)N_chord;         // timestep
    T scale = dx*dx*dx / (dt*dt);     // dx^3/dt^2 force scaling

    // Raw momentum sum -> physical force (factor 2 = momentum exchange)
    T Fx_phys = 2.0 * Fx * scale;
    T Fy_phys = 2.0 * Fy * scale;

    T q_inf = 0.5 * 1.0 * 1.0;        // rho=1, u_ref=1 (physical units)
    T area  = 1.0;                     // chord = 1

    T lift = -Fx_phys * std::sin(AoA_rad) + Fy_phys * std::cos(AoA_rad);
    T drag =  Fx_phys * std::cos(AoA_rad) + Fy_phys * std::sin(AoA_rad);

    Cl = lift / (q_inf * area);
    Cd = drag / (q_inf * area);
}

void writeVTK(MultiGridLattice2D<T, DESCRIPTOR>& lattice, plint iter, Box2D const& refineBox) {
    for (plint iLevel = 0; iLevel < lattice.getNumLevels(); ++iLevel) {
        MultiBlockLattice2D<T, DESCRIPTOR>& c = lattice.getComponent(iLevel);

        // Level 0: write the full domain. Level >0: write only the real
        // refined patch (its reported bounding box is the WHOLE domain
        // scaled by 2^level, but only refineBox*2^level actually has data).
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
//Simulation Loop
//===================================================================================================
void runSimulation(MultiGridLattice2D<T, DESCRIPTOR>& lattice, Box2D const& refineBox,
                    plint interpLevel, Array<plint,2> forceIds)
{
    using namespace param;
    //Open force output file
    std::ofstream forceFile("forces.txt");
    forceFile << "iter,avgEnergy\n";

    //Convergence Monitor
util::ValueTracer<T> convergence(param::U_lb, (T)param::N_chord, param::convTol);

    pcout << "Starting simulation: " << maxIter << " max iterations\n";

    for (plint iT = 0; iT <= maxIter; ++iT) {

    // Cheap: log forces every iteration (or every few) for proper time-resolution
    if (iT % forceLogIter == 0) {
        T Fx = lattice.getComponent(interpLevel).getInternalStatistics().getSum(forceIds[0]);
        T Fy = lattice.getComponent(interpLevel).getInternalStatistics().getSum(forceIds[1]);
        T Cl, Cd;
        computeCoefficients(Fx, Fy, Cl, Cd);
        forceFile << iT << "," << Cl << "," << Cd << "\n";
    }

    // Expensive: energy/VTK/convergence check at coarser cadence
    if (iT % outIter == 0) {
    T avgEnergy = computeAverageEnergy(lattice.getComponent(0));
    pcout << "iter=" << iT << " E=" << avgEnergy << "\n";
    convergence.takeValue(avgEnergy, true);
    if (convergence.hasConverged()) {
        pcout << "Converged at iter=" << iT << "\n";
        writeVTK(lattice, iT, refineBox);   // capture final state if converged early
        break;
        }
    }

    if (iT == 0 || iT == maxIter/3 || iT == 2*maxIter/3 || iT == maxIter) {
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

    //Startup Summary
    pcout << "============================================\n";
    pcout << "   NACA 0012 LBM Simulation - Palabos 2D   \n";
    pcout << "============================================\n";
    pcout << "Re        = " << param::Re       << "\n";
    pcout << "AoA       = " << param::AoA_deg  << " deg\n";
    pcout << "Domain    = " << param::Lx << " x " << param::Ly << " cells\n";
    pcout << "Chord     = " << param::N_chord  << " cells\n";
    pcout << "tau       = " << param::tau      << "\n";
    pcout << "omega     = " << param::omega    << "\n";
    pcout << "U_lb      = " << param::U_lb     << "\n";
    pcout << "Max iters = " << param::maxIter  << "\n";
    pcout << "Collision  = "
      << (param::collisionModel == param::BGK ? "BGK" : "MRT") << "\n";
    pcout << "============================================\n";
    

    //Write airfoil geometry CSV
    writeAirfoilGeometry("airfoil_geometry.csv");
    pcout << "Airfoil geometry written to airfoil_geometry.csv\n";

    //Build flag matrix and lattice
plint numLevel = 1, overlapWidth = 1, behaviorLevel = 0;
    MultiGridManagement2D management(param::Lx, param::Ly, numLevel, overlapWidth, behaviorLevel);
    Box2D refineBox(param::x_foil - param::N_chord/2, param::x_foil + 2*param::N_chord,
                    param::y_foil - param::N_chord,   param::y_foil + param::N_chord);
    // NOTE: management.refine(0, refineBox) removed — no refinement for this single-level pass

    plint xTiles = global::mpi().getSize(), yTiles = 1, interpLevel = numLevel - 1;   // interpLevel = 0
    ParallellizeBySquares2D* parallelizer = new ParallellizeBySquares2D(
        management.getBulks(), management.getBoundingBox(interpLevel), xTiles, yTiles);
    management.parallelize(parallelizer);

    MultiGridLattice2D<T, DESCRIPTOR> lattice =
        *MultiGridGenerator2D<T, DESCRIPTOR>::createRefinedLatticeCubicInterpolationFiltering(
            management, makeBulkDynamics(), behaviorLevel);
    pcout << "Lattice built: " << param::Lx * param::Ly << " cells\n";

    Box2D bb0 = lattice.getComponent(0).getBoundingBox();
    pcout << "Level 0 bbox: x[" << bb0.x0 << "," << bb0.x1 << "] y[" << bb0.y0 << "," << bb0.y1 << "]\n";
    // Level 1 bbox print removed — no level 1 exists in this configuration

    // airfoilDomain: global (level-0) coordinates. interpLevel == 0 here, so no rescaling needed.
    Box2D airfoilDomain(param::x_foil - 2, param::x_foil + param::N_chord + 2,
                         param::y_foil - param::N_chord/4, param::y_foil + param::N_chord/4);
    
    Array<plint,2> forceIds;
    forceIds[0] = lattice.getComponent(interpLevel).internalStatSubscription().subscribeSum();
    forceIds[1] = lattice.getComponent(interpLevel).internalStatSubscription().subscribeSum();


// Stamp airfoil geometry on every grid level
for (plint iLevel = 0; iLevel < numLevel; ++iLevel) {
    stampAirfoilOnLevel(lattice, iLevel, forceIds);
}

initializeMomentumExchange(lattice.getComponent(interpLevel), airfoilDomain);

    //Boundary Conditions
    setBoundaryConditions(lattice);
    pcout << "Boundary conditions applied\n";

    //Initialise distribution functions
    initializeDomain(lattice);
    pcout << "Domain initialised at equilibrium\n";

    initializeDomain(lattice);
    pcout << "Domain initialised at equilibrium\n";
    pcout << "t=0 check: avgEnergy = " << computeAverageEnergy(lattice.getComponent(0)) << "\n";

    //run
    runSimulation(lattice, refineBox, interpLevel, forceIds);
    
    return 0;
}