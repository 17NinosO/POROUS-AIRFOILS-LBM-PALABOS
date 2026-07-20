// porous_channel_mg.cpp
#include "palabos2D.h"
#include "palabos2D.hh"

#include <csignal>
#include <iostream>
#include <cmath>

#include "config.hpp"
#include "domain.hpp"
#include "pore_mask_domain.hpp"

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::MRTD2Q9Descriptor

namespace param {
    const T targetRe = 50.0;
    const T U_lb = 0.02;
    const plint maxIter = 200000;
    const plint checkIter = 500;
    const T convergenceThreshold = 1e-4;
}

namespace {
    volatile std::sig_atomic_t g_stopRequested = 0;
    void handleSigint(int) { g_stopRequested = 1; }
}

class SetPoreSolidOnLevel : public BoxProcessingFunctional2D_S<int> {
public:
    SetPoreSolidOnLevel(PoreMaskDomain domain_, plint blockX0_, plint blockY0_)
        : domain(std::move(domain_)), blockX0(blockX0_), blockY0(blockY0_) {}
    virtual void process(Box2D box, ScalarField2D<int>& flags) {
        Dot2D offset = flags.getLocation();
        for (plint ix = box.x0; ix <= box.x1; ++ix) {
            for (plint iy = box.y0; iy <= box.y1; ++iy) {
                plint gx = ix + offset.x, gy = iy + offset.y;
                plint bx = gx - blockX0, by = gy - blockY0;
                if (bx < 0 || bx >= domain.nx() || by < 0 || by >= domain.ny()) continue;
                if (domain.isSolid(static_cast<int>(bx), static_cast<int>(by)))
                    flags.get(ix, iy) = 1;
            }
        }
    }
    virtual SetPoreSolidOnLevel* clone() const { return new SetPoreSolidOnLevel(*this); }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
    }
private:
    PoreMaskDomain domain;
    plint blockX0, blockY0;
};

void stampPoreOnLevel(MultiGridLattice2D<T, DESCRIPTOR>& lattice, plint iLevel,
                       std::vector<PoreMaskDomain>& levelDomains, plint blockX0, plint blockY0) {
    MultiBlockLattice2D<T, DESCRIPTOR>& comp = lattice.getComponent(iLevel);
    plint scale = 1 << iLevel;
    MultiScalarField2D<int> flags(comp);
    Box2D bb = comp.getBoundingBox();
    applyProcessingFunctional(
        new SetPoreSolidOnLevel(levelDomains[iLevel], blockX0 * scale, blockY0 * scale), bb, flags);
    defineDynamics(comp, flags, new BounceBack<T, DESCRIPTOR>(), 1);
    pcout << "Pore stamped on level " << iLevel << "\n";
}

// Unmodified from porous_channel.cpp -- applied to a single component
// (the coarse level), not the multigrid lattice as a whole.
void setBoundaryConditions(MultiBlockLattice2D<T, DESCRIPTOR>& lattice, plint Lx, plint Ly) {
    Array<T, 2> u_inf(param::U_lb, 0.0);
    OnLatticeBoundaryCondition2D<T, DESCRIPTOR>* bc = createLocalBoundaryCondition2D<T, DESCRIPTOR>();
    Box2D inlet(0, 0, 1, Ly - 2);
    Box2D outlet(Lx - 1, Lx - 1, 1, Ly - 2);
    Box2D bottomWall(0, Lx - 1, 0, 0);
    Box2D topWall(0, Lx - 1, Ly - 1, Ly - 1);
    bc->setVelocityConditionOnBlockBoundaries(lattice, inlet);
    bc->setVelocityConditionOnBlockBoundaries(lattice, bottomWall);
    bc->setVelocityConditionOnBlockBoundaries(lattice, topWall);
    bc->setVelocityConditionOnBlockBoundaries(lattice, outlet, boundary::outflow);
    setBoundaryVelocity(lattice, lattice.getBoundingBox(), u_inf);
    delete bc;
}

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);
    std::signal(SIGINT, handleSigint);

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
    config.resolution = 1.0;

    auto mask = pore_geometry::buildDomainMask(config);
    PoreMaskDomain poreDomain(mask);

    std::vector<PoreMaskDomain> levelDomains;
    for (plint iLevel = 0; iLevel < 4; ++iLevel) {
        pore_geometry::GeometryConfig levelConfig = config;
        levelConfig.resolution = config.resolution / static_cast<double>(1 << iLevel);
        levelDomains.emplace_back(pore_geometry::buildDomainMask(levelConfig));
    }

    plint blockNx = poreDomain.nx();
    plint blockNy = poreDomain.ny();
    plint Lx = 10 * blockNx;
    plint Ly = 10 * blockNy;
    plint blockX0 = (Lx - blockNx) / 3;
    plint blockY0 = (Ly - blockNy) / 2;

    T L_lb = static_cast<T>(blockNy);
    T nu_lb = param::U_lb * L_lb / param::targetRe;
    T tau = 3.0 * nu_lb + 0.5;
    T omega = 1.0 / tau;
    pcout << "Re=" << param::targetRe << " tau=" << tau << "\n";
    pcout << "Channel = " << Lx << " x " << Ly << " (block " << blockNx << " x " << blockNy << ")\n";

    plint numLevel = 4, overlapWidth = 1, behaviorLevel = 0;
    MultiGridManagement2D management(Lx, Ly, numLevel, overlapWidth, behaviorLevel);
    plint margin = 20;
    Box2D refineBox(blockX0 - margin, blockX0 + blockNx - 1 + margin,
                     blockY0 - margin, blockY0 + blockNy - 1 + margin);
    management.refine(0, refineBox);
    management.refine(1, refineBox.multiply(2));
    management.refine(2, refineBox.multiply(4));

    plint xTiles = global::mpi().getSize(), yTiles = 1, interpLevel = numLevel - 1;
    ParallellizeBySquares2D* parallelizer = new ParallellizeBySquares2D(
        management.getBulks(), management.getBoundingBox(interpLevel), xTiles, yTiles);
    management.parallelize(parallelizer);

    MultiGridLattice2D<T, DESCRIPTOR> lattice =
        *MultiGridGenerator2D<T, DESCRIPTOR>::createRefinedLatticeCubicInterpolationFiltering(
            management, new MRTdynamics<T, DESCRIPTOR>(omega), behaviorLevel);

    for (plint iLevel = 0; iLevel < lattice.getNumLevels(); ++iLevel)
        stampPoreOnLevel(lattice, iLevel, levelDomains, blockX0, blockY0);

    setBoundaryConditions(lattice.getComponent(0), Lx, Ly);

    Array<T, 2> u0(param::U_lb, 0.0);
    for (plint iLevel = 0; iLevel < lattice.getNumLevels(); ++iLevel) {
        MultiBlockLattice2D<T, DESCRIPTOR>& comp = lattice.getComponent(iLevel);
        initializeAtEquilibrium(comp, comp.getBoundingBox(), (T)1.0, u0);
    }
    lattice.initialize();

    // Convergence checked on the finest level (3), where the actual
    // pore-scale flow lives, restricted to the block's own footprint
    // scaled into level-3 coordinates.
    plint finestScale = 1 << 3;
    Box2D blockRegionFinest(blockX0 * finestScale, (blockX0 + blockNx - 1) * finestScale,
                             blockY0 * finestScale, (blockY0 + blockNy - 1) * finestScale);

    plb_ofstream dataFile("porous_channel_mg_data.csv");
    dataFile << "iter,avgEnergy\n";

    T previousEnergy = -1.0;
    for (plint iT = 0; iT <= param::maxIter && !g_stopRequested; ++iT) {
        if (iT % param::checkIter == 0) {
            T energy = computeAverageEnergy(lattice.getComponent(3), blockRegionFinest);
            if (!std::isfinite(energy)) {
                pcout << "iter=" << iT << " avgEnergy=" << energy << "  <- DIVERGED, stopping\n";
                break;
            }
            pcout << "iter=" << iT << " avgEnergy=" << energy;
            dataFile << iT << "," << energy << "\n";

            if (previousEnergy > 0.0) {
                T relChange = std::abs(energy - previousEnergy) / previousEnergy;
                pcout << " relChange=" << relChange;
                if (relChange < param::convergenceThreshold) {
                    pcout << "  <- CONVERGED\n";
                    break;
                }
            }
            pcout << "\n";
            previousEnergy = energy;

            Box2D refineBoxForWrite(blockX0, blockX0 + blockNx - 1, blockY0, blockY0 + blockNy - 1);
            for (plint iLevel = 0; iLevel < lattice.getNumLevels(); ++iLevel) {
                MultiBlockLattice2D<T, DESCRIPTOR>& comp = lattice.getComponent(iLevel);
                plint scale = 1 << iLevel;
                Box2D writeBox = comp.getBoundingBox();
                if (iLevel > 0)
                    writeBox = Box2D(refineBoxForWrite.x0 * scale, refineBoxForWrite.x1 * scale,
                                      refineBoxForWrite.y0 * scale, refineBoxForWrite.y1 * scale);
                VtkImageOutput2D<T> vtkOut(
                    createFileName("porous_mg_level" + std::to_string(iLevel), iT, 6),
                    1.0 / static_cast<double>(scale));
                vtkOut.writeData<2, float>(*computeVelocity(comp, writeBox), "velocity", 1.0);
            }
        }
        lattice.collideAndStream();
    }

    pcout << "Run complete.\n";
    return 0;
}