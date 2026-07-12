// porous_channel.cpp
#include "palabos2D.h"
#include "palabos2D.hh"

#include <iostream>

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::D2Q9Descriptor

namespace param {
    const plint Lx = 200;   // placeholder domain size -- Step 2c will size
    const plint Ly = 100;   // this from the actual pore_geometry mask instead
    const T tau = 0.8;      // placeholder -- Step 3 derives this from Re
    const T omega = 1.0 / tau;
    const T U_lb = 0.02;    // small, safe placeholder inlet-ish velocity
    const plint maxIter = 500;
    const plint outIter = 100;
}

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