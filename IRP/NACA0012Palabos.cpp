#include "palabos2D.h"  //Declarations:class interfaces, types
#include "palabos2D.hh" //Definitions: template implementations

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>

#include "complexDynamics/mrtDynamics.h"
#include "complexDynamics/mrtDynamics.hh"

using namespace plb;    //Palabos namespace
using namespace std;

// Palabos is a template library (e.g. MultiBlockLattice2D<T, DESCRIPTOR>).
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
    const T AoA_deg = 5.0; //Angle of attack [degrees]
    const T AoA_rad = AoA_deg * M_PI / 180.0; //Angle of attack [radians]

    //Lattice Resolution
    //N_chord: number of lattice cells spanning the chord
    //Higher = more accurate but more computationally expensive
    //100 is a good starting point; but use 50 for a quick test run.
    const plint N_chord = 200;

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
    const plint Lx = 25 * N_chord; // Domain Width
    const plint Ly = 15 * N_chord; // Domain Height 
    const plint x_foil = 8 * N_chord; // LE x position from the inlet.
    const plint y_foil = Ly / 2; // LE y centred vertically.

    //Simulation Control Parameters
    const plint maxIter = 200000; //Maximum number of iterations.
    const plint outIter = 1000; //Output and Logging Interval
    const T convTol = 1e-6; //Velocity convergence tolerance.

    // Collision model
    enum CollisionModel { BGK, MRT };
    const CollisionModel collisionModel = MRT;   //change BGK/MRT

    //Porosity Parameters
    const bool porous = true;
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

bool isInsideAirfoil(plint ix, plint iy) {
    using namespace param;

    //Map to normalised chord frame
    T x_norm = static_cast<T>(ix - x_foil) / static_cast<T>(N_chord);
    T y_norm = static_cast<T>(iy - y_foil) / static_cast<T>(N_chord);

    //Outside chord extent - definitely fluid
    if (x_norm < 0.0 || x_norm > 1.0) return false;

    //Compute half thickness at this chord position
    T yt = naca0012Thickness(x_norm);

    bool insideBody = (std::fabs(y_norm) < yt);

    return insideBody && !isInsidePore(ix, iy);
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
//Creating the MultiBlockLattice and stamping the correct collision dynamics onto each cell type.

namespace cellFlag {
    const int FLUID = 0;    //Regular fluid cell
    const int SOLID = 1;    //Airfoil interior
    const int INLET = 2;    //Left boundary
    const int OUTLET = 3;   //Right boundary
    const int WALL = 4;     //Top and bottom boundaries
}

//Palabos decomposes the domain into blocks for parallel processing.
//When you want to loop over cells you can't write a flat nested for loop over global indices.
//Each parallel process only owns a portoin of the domain and the loop coordinates are local to that block.

//Palabos solves this with processing functionals
//These are objects that receive a block plus its local domain.
//They use a getLocation() to convert local into global coordinates.

//Stamping of solid flags onto all cells inside the airfoil geometry.
//Inherits from BoxProcessingFunctional2D_S<int> - the <int> matches
//the scaler field type (MultiScalarField2D<int>).

class SetAirfoilSolid : public BoxProcessingFunctional2D_S<int> {
public:
    void process(Box2D domain, ScalarField2D<int>& flags) {

        //getLocation() returns this blocks offset in the global coordinates.
        //without this, ix/iy are local - wrong geometry placement.
        Dot2D offset = flags.getLocation();

        for (plint ix = domain.x0; ix <= domain.x1; ++ix) {
            for (plint iy = domain.y0; iy<= domain.y1; ++iy) {

                //Convert to global lattice coordinates
                plint gx = ix + offset.x;
                plint gy = iy + offset.y;

                if (isInsideAirfoil(gx, gy)) {
                    flags.get(ix, iy) = cellFlag::SOLID;
                }
            }
        }
    }

    //Palabos requires these two boilerplate mathods on every functional:
    //clone() lets Palabos copy the functional for each block it processes.
    SetAirfoilSolid* clone() const {
        return new SetAirfoilSolid(*this);
    }

    //Tell Palabos that this is written to the scalar field.
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
    }
};

//Using setToConstant for the simple rectangular regions and the custom functional for the airfoil.
//Populate a pre-allocated MultiScalarField2D<int> with cell type flags

void buildGeometryFlags(MultiScalarField2D<int>& flags) {
    using namespace param;

    //Have every cell start as a fluid
    setToConstant(flags, flags.getBoundingBox(), cellFlag::FLUID);

    //Boundary Strips - Box2D(x0, x1, y0, y1) defines a rectangle region
    setToConstant(flags, Box2D(0, 0, 0, Ly-1), cellFlag::INLET); //Left strip
    setToConstant(flags, Box2D(Lx-1, Lx-1, 0, Ly-1), cellFlag::OUTLET); //Right strip
    setToConstant(flags, Box2D(0, Lx-1, 0, 0), cellFlag::WALL); //Bottom strip
    setToConstant(flags, Box2D(0, Lx-1, Ly-1, Ly-1), cellFlag::WALL); //Top strip

    //Stamp the airfoil solid cells using the custom functional
    //applyProcessingFunctional dispatches across the blocks
    applyProcessingFunctional(
        new SetAirfoilSolid(),
        flags.getBoundingBox(),
        flags
    );
}

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
//Create the MultiBlockLattice2D, and assign dynamics from the flag matrix.
//Returns a raw pointer

MultiBlockLattice2D<T, DESCRIPTOR>* buildLattice(
    MultiScalarField2D<int>& flags,
    Array<plint,2>& forceIds)                    // <-- NEW parameter
{
    using namespace param;

    MultiBlockLattice2D<T, DESCRIPTOR>* lattice =
        new MultiBlockLattice2D<T, DESCRIPTOR>(
            Lx, Ly, makeBulkDynamics());

    lattice->periodicity().toggleAll(false);

    //Subscribe x and y force components to internal statistics
    forceIds[0] = lattice->internalStatSubscription().subscribeSum();
    forceIds[1] = lattice->internalStatSubscription().subscribeSum();

    //MomentumExchangeBounceBack instead of plain BounceBack
    defineDynamics(*lattice, flags,
                   new MomentumExchangeBounceBack<T, DESCRIPTOR>(forceIds),
                   cellFlag::SOLID);

    // NEW: required to initialise the momentum-exchange machinery
    initializeMomentumExchange(*lattice, lattice->getBoundingBox());

    return lattice;
}

//===================================================================================================
//Initialisation
//===================================================================================================

void initializeDomain(MultiBlockLattice2D<T, DESCRIPTOR>& lattice) {
    using namespace param;

    //Initial Density
    //In LBM lattice units, density is normalised: rho = 1
    //Pressure is related by: rho = rho * cs_sq, so p = rho/3
    const T rho_0 = 1.0;

    //Free stream velocity
    //This is the AoA implementation.

    Array<T,2> u_inf(
        U_lb * std::cos(AoA_rad),
        U_lb * std::sin(AoA_rad)
    );

    //Set the distribution function f_i = f_ieq(rho_0, u_inf) in the domain
    //Palabos computes the full equilibrium formula internally.
    //We pass the bounding box so it applies to the entire lattice.

    initializeAtEquilibrium(
        lattice,
        lattice.getBoundingBox(),
        rho_0,
        u_inf
    );

    lattice.initialize();
}

//===================================================================================================
//Boundary Conditions
//===================================================================================================
void setBoundaryConditions(MultiBlockLattice2D<T, DESCRIPTOR>& lattice) {
    using namespace param;

    //Free-stream velocity vector where AoA is encoded as direction
    Array<T, 2> u_inf(
        U_lb * std::cos(AoA_rad),
        U_lb * std::sin(AoA_rad)
    );

    //Create the Zou-He boundary condition object
    OnLatticeBoundaryCondition2D<T, DESCRIPTOR>* bc =
    createLocalBoundaryCondition2D<T, DESCRIPTOR>();

    //Inlet Left wall
    //Skip corner cells
    bc->addVelocityBoundary0N(Box2D(0, 0, 1, Ly-2), lattice);
    setBoundaryVelocity(lattice, Box2D(0, 0, 1, Ly-2), u_inf);

    //Outlet Right wall
    //Density as 1
    //Velocity not prescribed to let flow exit freely
    bc->addPressureBoundary0P(Box2D(Lx-1, Lx-1, 1, Ly-2), lattice);
    setBoundaryDensity(lattice, Box2D(Lx-1, Lx-1, 1, Ly-2), (T)1.0);

    //Top wall
    //Free stream velocity condition
    //Domain is 15c tall so is reasonable
    bc->addVelocityBoundary1P(Box2D(1, Lx-2, Ly-1, Ly-1), lattice);
    setBoundaryVelocity(lattice, Box2D(1, Lx-2, Ly-1, Ly-1), u_inf);

    //Bottom wall
    bc->addVelocityBoundary1N(Box2D(1, Lx-2, 0, 0), lattice);
    setBoundaryVelocity(lattice, Box2D(1, Lx-2, 0, 0), u_inf);

    //Corner Walls
    defineDynamics(lattice, Box2D(0,    0,    0,    0   ), new BounceBack<T,DESCRIPTOR>());  // bottom-left
    defineDynamics(lattice, Box2D(0,    0,    Ly-1, Ly-1), new BounceBack<T,DESCRIPTOR>());  // top-left
    defineDynamics(lattice, Box2D(Lx-1, Lx-1, 0,    0   ), new BounceBack<T,DESCRIPTOR>());  // bottom-right
    defineDynamics(lattice, Box2D(Lx-1, Lx-1, Ly-1, Ly-1), new BounceBack<T,DESCRIPTOR>());  // top-right

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

void writeVTK(MultiBlockLattice2D<T, DESCRIPTOR>& lattice, plint iter) {
    VtkImageOutput2D<T> vtkOut(createFileName("vtk", iter, 6), 1.0);
    vtkOut.writeData<float>(
        *computeVelocityNorm(lattice), "velocityNorm", 1.0);
    vtkOut.writeData<2, float>(
        *computeVelocity(lattice), "velocity", 1.0);
    vtkOut.writeData<float>(
        *computeDensity(lattice), "density", 1.0);
}

//===================================================================================================
//Simulation Loop
//===================================================================================================
void runSimulation(
    MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
    MultiScalarField2D<int>& flags,
    Array<plint,2>& forceIds)
{
    using namespace param;
    //Open force output file
    std::ofstream forceFile("forces.txt");
    forceFile << "iter,Cl,Cd,avgEnergy\n";

    //Convergence Monitor
util::ValueTracer<T> convergence(param::U_lb, (T)param::N_chord, param::convTol);

    pcout << "Starting simulation: " << maxIter << " max iterations\n";

    for (plint iT = 0; iT <= maxIter; ++iT) {
        //Output block
        if (iT % outIter == 0) {
            T avgEnergy = computeAverageEnergy(lattice);

            if (iT > 5000) {
                T Fx = lattice.getInternalStatistics().getSum(forceIds[0]);
                T Fy = lattice.getInternalStatistics().getSum(forceIds[1]);
                T Cl, Cd;
                computeCoefficients(Fx, Fy, Cl, Cd);

                forceFile << iT << "," << Cl << "," << Cd << ","
                          << avgEnergy << "\n";

                pcout << "iter=" << iT
                    << " Cl=" << Cl
                    << " Cd=" << Cd
                    << " E=" << avgEnergy << "\n";
            } else {
                pcout << "iter=" << iT
                    << " [transient] E=" << avgEnergy << "\n";
            }

            //Write VTK snapshot
            writeVTK(lattice, iT);

            //Feed convergence monitor
            convergence.takeValue(avgEnergy, true);
            if (convergence.hasConverged()) {
                pcout << "Converged at iter=" << iT << "\n";
                break;
            }
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
MultiScalarField2D<int> flags(param::Lx, param::Ly);
    buildGeometryFlags(flags);
    Array<plint,2> forceIds;
    MultiBlockLattice2D<T, DESCRIPTOR>* lattice = buildLattice(flags, forceIds);
    pcout << "Lattice built: " << param::Lx * param::Ly << " cells\n";

    //Boundary Conditions
    setBoundaryConditions(*lattice);
    pcout << "Boundary conditions applied\n";

    //Initialise distribution functions
    initializeDomain(*lattice);
    pcout << "Domain initialised at equilibrium\n";

    //run
runSimulation(*lattice, flags, forceIds);

    delete lattice;
    return 0;
}