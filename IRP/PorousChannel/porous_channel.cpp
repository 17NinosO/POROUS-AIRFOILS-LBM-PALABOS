// porous_channel.cpp
#include "palabos2D.h"
#include "palabos2D.hh"

#include <iostream>

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::D2Q9Descriptor

namespace param {
    const plint Lx = 200;
    const plint Ly = 100;
    const T tau = 0.8;
    const T omega = 1.0 / tau;
    const T U_lb = 0.02;
    const plint maxIter = 500;
    const plint outIter = 100;
}

// -------------------- NEW: goes here, above main() --------------------
void setBoundaryConditions(MultiBlockLattice2D<T, DESCRIPTOR>& lattice) {
    using namespace param;
    Array<T, 2> u_inf(U_lb, 0.0);

    OnLatticeBoundaryCondition2D<T, DESCRIPTOR>* bc =
        createLocalBoundaryCondition2D<T, DESCRIPTOR>();

    bc->addVelocityBoundary0N(Box2D(0, 0, 1, Ly - 2), lattice);
    setBoundaryVelocity(lattice, Box2D(0, 0, 1, Ly - 2), u_inf);

    bc->addPressureBoundary0P(Box2D(Lx - 1, Lx - 1, 1, Ly - 2), lattice);
    setBoundaryDensity(lattice, Box2D(Lx - 1, Lx - 1, 1, Ly - 2), (T)1.0);

    defineDynamics(lattice, Box2D(1, Lx - 2, 0, 0), new BounceBack<T, DESCRIPTOR>());
    defineDynamics(lattice, Box2D(1, Lx - 2, Ly - 1, Ly - 1), new BounceBack<T, DESCRIPTOR>());

    defineDynamics(lattice, Box2D(0, 0, 0, 0), new BounceBack<T, DESCRIPTOR>());
    defineDynamics(lattice, Box2D(0, 0, Ly - 1, Ly - 1), new BounceBack<T, DESCRIPTOR>());
    defineDynamics(lattice, Box2D(Lx - 1, Lx - 1, 0, 0), new BounceBack<T, DESCRIPTOR>());
    defineDynamics(lattice, Box2D(Lx - 1, Lx - 1, Ly - 1, Ly - 1), new BounceBack<T, DESCRIPTOR>());

    delete bc;
}
// ------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);

    pcout << "==========================================\n";
    pcout << "  Porous channel -- build/link smoke test  \n";
    pcout << "==========================================\n";
    pcout << "Domain = " << param::Lx << " x " << param::Ly << "\n";
    pcout << "tau    = " << param::tau << "  (placeholder)\n";

    MultiBlockLattice2D<T, DESCRIPTOR> lattice(
        param::Lx, param::Ly,
        new BGKdynamics<T, DESCRIPTOR>(param::omega));

    setBoundaryConditions(lattice);   // <-- NEW: this one line, right here

    Array<T, 2> u0(param::U_lb, 0.0);
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), (T)1.0, u0);
    lattice.initialize();

    pcout << "t=0 avgEnergy = " << computeAverageEnergy(lattice) << "\n";

    for (plint iT = 0; iT <= param::maxIter; ++iT) {
        if (iT % param::outIter == 0) {
            pcout << "iter=" << iT << " avgEnergy=" << computeAverageEnergy(lattice) << "\n";
            VtkImageOutput2D<T> vtkOut(createFileName("smoke_test", iT, 6), 1.0);
            vtkOut.writeData<2, float>(*computeVelocity(lattice), "velocity", 1.0);
        }
        lattice.collideAndStream();
    }

    pcout << "Smoke test complete -- build chain confirmed working.\n";
    return 0;
}