// porous_channel.cpp
#include "palabos2D.h"
#include "palabos2D.hh"

#include <iostream>

#include "config.hpp"
#include "domain.hpp"
#include "pore_mask_domain.hpp"

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::MRTD2Q9Descriptor

namespace param {
    const T targetRe = 500.0;
    const T U_lb = 0.02;
    const plint maxIter = 200000;
    const plint checkIter = 200;
    const T convergenceThreshold = 1e-4;
}

// Mirrors Gabriel's MaskShapeDomain2D: a DomainFunctional2D predicate,
// which is what MomentumExchangeBounceBack + initializeMomentumExchange
// require -- not a pre-built flags field.
class PoreBlockDomain : public DomainFunctional2D {
public:
    PoreBlockDomain(PoreMaskDomain domain_, plint blockX0_, plint blockY0_, bool solidOnly_ = false)
        : domain(std::move(domain_)), blockX0(blockX0_), blockY0(blockY0_), solidOnly(solidOnly_) {}

    virtual bool operator()(plint iX, plint iY) const {
        plint bx = iX - blockX0;
        plint by = iY - blockY0;
        if (bx < 0 || bx >= domain.nx() || by < 0 || by >= domain.ny())
            return false;  // outside the block -> freestream, not solid
        if (solidOnly) return true;  // validation mode: plain solid block, no pores
        return domain.isSolid(static_cast<int>(bx), static_cast<int>(by));
    }
    
private:
    PoreMaskDomain domain;
    plint blockX0, blockY0;
    bool solidOnly;

    virtual PoreBlockDomain* clone() const { return new PoreBlockDomain(*this); }

};

void stampPoreGeometry(MultiBlockLattice2D<T, DESCRIPTOR>& lattice, const PoreMaskDomain& domain,
                        plint blockX0, plint blockY0, Array<plint, 2> forceIds, bool solidOnly) {
    Box2D bb = lattice.getBoundingBox();
    defineDynamics(
        lattice, bb, new PoreBlockDomain(domain, blockX0, blockY0, solidOnly),
        new MomentumExchangeBounceBack<T, DESCRIPTOR>(forceIds));
    initializeMomentumExchange(
        lattice, bb, new PoreBlockDomain(domain, blockX0, blockY0, solidOnly));
    pcout << "Pore block stamped (momentum exchange enabled) at (" << blockX0 << ", " << blockY0 << ")\n";
}

void setBoundaryConditions(MultiBlockLattice2D<T, DESCRIPTOR>& lattice, plint Lx, plint Ly) {
    using namespace param;
    Array<T, 2> u_inf(U_lb, 0.0);

    OnLatticeBoundaryCondition2D<T, DESCRIPTOR>* bc =
        createLocalBoundaryCondition2D<T, DESCRIPTOR>();

    Box2D inlet(0, 0, 1, Ly - 2);
    Box2D outlet(Lx - 1, Lx - 1, 1, Ly - 2);
    Box2D bottomWall(0, Lx - 1, 0, 0);
    Box2D topWall(0, Lx - 1, Ly - 1, Ly - 1);

    // Top/bottom treated as inlets, matching Gabriel's proven approach --
    // not classic no-slip bounce-back walls.
    bc->setVelocityConditionOnBlockBoundaries(lattice, inlet);
    bc->setVelocityConditionOnBlockBoundaries(lattice, bottomWall);
    bc->setVelocityConditionOnBlockBoundaries(lattice, topWall);
    bc->setVelocityConditionOnBlockBoundaries(lattice, outlet, boundary::outflow);

    setBoundaryVelocity(lattice, lattice.getBoundingBox(), u_inf);

    delete bc;
}

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);

    bool solidOnly = false;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--solid") solidOnly = true;
    }
    pcout << "Mode: " << (solidOnly ? "SOLID BLOCK (validation)" : "porous geometry") << "\n";

    pcout << "==========================================\n";
    pcout << "     Porous channel -- pore geometry run   \n";
    pcout << "==========================================\n";

    // Build the geometry FIRST -- its resolution decides the lattice size
    pore_geometry::GeometryConfig config;
    config.disorder = 0.0;
    config.r_mean = 10.0;
    config.r_std = 1.0;
    config.r_min = 8.0;
    config.r_max = 12.0;
    config.throat_width_mean = 3.0;
    config.throat_width_std = 0.5;
    config.throat_width_min = 1.0;
    config.boundary_roughness = 1.0;
    config.resolution = 1.0;   // coarser than default -- faster first run

    auto mask = pore_geometry::buildDomainMask(config);
    PoreMaskDomain poreDomain(mask);

    plint blockNx = poreDomain.nx();
    plint blockNy = poreDomain.ny();

    plint Lx = 10 * blockNx;
    plint Ly = 10 * blockNy;
    plint blockX0 = (Lx - blockNx) / 3;
    plint blockY0 = (Ly - blockNy) / 2;

    // Re -> tau: L_lb is the block's cross-stream height, the
    // characteristic length for the wake/bluff-body flow we're watching
    T L_lb = static_cast<T>(blockNy);
    T nu_lb = param::U_lb * L_lb / param::targetRe;
    T tau = 3.0 * nu_lb + 0.5;
    T omega = 1.0 / tau;
    pcout << "Re=" << param::targetRe << " L_lb=" << L_lb << " nu_lb=" << nu_lb << " tau=" << tau << "\n";

    Box2D blockRegion(blockX0, blockX0 + blockNx - 1, blockY0, blockY0 + blockNy - 1);

    pcout << "Channel = " << Lx << " x " << Ly << " (block " << blockNx << " x " << blockNy << ")\n";

    MultiBlockLattice2D<T, DESCRIPTOR> lattice(
        Lx, Ly, new MRTdynamics<T, DESCRIPTOR>(omega));

// Register two internal-statistics sums BEFORE stamping -- this is
    // what forceIds actually points at; MomentumExchangeBounceBack
    // accumulates into these during every collideAndStream().
    Array<plint, 2> forceIds;
    forceIds[0] = lattice.internalStatSubscription().subscribeSum();
    forceIds[1] = lattice.internalStatSubscription().subscribeSum();

    stampPoreGeometry(lattice, poreDomain, blockX0, blockY0, forceIds, solidOnly);
    setBoundaryConditions(lattice, Lx, Ly);

    // Same lattice -> physical scaling convention as the NACA0012 case,
    // with L_lb (block height) playing the role N_chord played there.
    T dx = 1.0 / L_lb;
    T dt = param::U_lb / L_lb;
    T scale = dx * dx * dx / (dt * dt);

    Array<T, 2> u0(param::U_lb, 0.0);
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), (T)1.0, u0);
    lattice.initialize();

    pcout << "t=0 avgEnergy = " << computeAverageEnergy(lattice, blockRegion) << "\n";

    plb_ofstream dataFile("porous_channel_data.csv");
    dataFile << "iter,avgEnergy,drag,lift,Cd,Cl\n";

    T previousEnergy = -1.0;
    for (plint iT = 0; iT <= param::maxIter; ++iT) {
        if (iT % param::checkIter == 0) {
            T energy = computeAverageEnergy(lattice, blockRegion);
            T Fx_phys = 2.0 * lattice.getInternalStatistics().getSum(forceIds[0]) * scale;
            T Fy_phys = 2.0 * lattice.getInternalStatistics().getSum(forceIds[1]) * scale;
            pcout << "iter=" << iT << " avgEnergy=" << energy
                  << " drag=" << Fx_phys << " lift=" << Fy_phys
                  << " Cd=" << Fx_phys / 0.5 << " Cl=" << Fy_phys / 0.5;
            dataFile << iT << "," << energy << "," << Fx_phys << "," << Fy_phys << ","
                      << Fx_phys / 0.5 << "," << Fy_phys / 0.5 << "\n";

            if (previousEnergy > 0.0) {
                T relChange = std::abs(energy - previousEnergy) / previousEnergy;
                pcout << " relChange=" << relChange;

                if (relChange < param::convergenceThreshold) {
                    pcout << "  <- CONVERGED\n";
                    VtkImageOutput2D<T> vtkOut(createFileName("porous_channel_final", iT, 6), 1.0);
                    vtkOut.writeData<2, float>(*computeVelocity(lattice), "velocity", 1.0);
                    break;
                }
            }
            pcout << "\n";
            previousEnergy = energy;

            VtkImageOutput2D<T> vtkOut(createFileName("porous_channel", iT, 6), 1.0);
            vtkOut.writeData<2, float>(*computeVelocity(lattice), "velocity", 1.0);
        }
        lattice.collideAndStream();
    }

    pcout << "Run complete.\n";
    return 0;
}