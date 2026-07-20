// multigrid_test.cpp
//
// Verifies MultiGridManagement2D construction, per-level pore mask
// regeneration, and per-level stamping together. Confirmed working:
// a full run through all 4 levels with correct stamping and VTK output.
//
// KNOWN OPEN ISSUE: an intermittent, nondeterministic crash (SIGFPE
// inside PoreShape::radiusAt, or a related heap issue) has been
// observed on some runs despite identical input/seed -- most likely
// an uninitialized-memory read somewhere in pore_geometry, not yet
// root-caused. If this recurs, don't assume it's fixed just because
// it ran clean once.

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

    MultiScalarField2D<int> flags(comp);  // shares comp's real block structure
    Box2D bb = comp.getBoundingBox();
    applyProcessingFunctional(
        new SetPoreSolidOnLevel(levelDomains[iLevel], blockX0 * scale, blockY0 * scale), bb, flags);
    defineDynamics(comp, flags, new BounceBack<T, DESCRIPTOR>(), 1);
    pcout << "Pore stamped on level " << iLevel << "\n";
}

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);

    // Real config first -- everything else depends on these values
    // being set BEFORE the block dimensions or per-level masks are
    // computed from them.
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

    // One DomainMask per level, same seed, each at that level's native
    // resolution -- built from the SAME fully-configured `config`, not
    // regenerated from defaults.
    std::vector<PoreMaskDomain> levelDomains;
    for (plint iLevel = 0; iLevel < 4; ++iLevel) {
        pore_geometry::GeometryConfig levelConfig = config;
        levelConfig.resolution = config.resolution / static_cast<double>(1 << iLevel);
        levelDomains.emplace_back(pore_geometry::buildDomainMask(levelConfig));
        pcout << "level " << iLevel << " mask: " << levelDomains.back().nx()
              << " x " << levelDomains.back().ny() << "\n";
    }

    plint blockNx = poreDomain.nx();
    plint blockNy = poreDomain.ny();
    plint Lx = 10 * blockNx;
    plint Ly = 10 * blockNy;
    plint blockX0 = (Lx - blockNx) / 3;
    plint blockY0 = (Ly - blockNy) / 2;

    pcout << "Coarse channel: " << Lx << " x " << Ly << "\n";
    pcout << "Block (level 0 coords): x0=" << blockX0 << " y0=" << blockY0
          << " nx=" << blockNx << " ny=" << blockNy << "\n";

    // --- multigrid construction ---
    plint numLevel = 4, overlapWidth = 1, behaviorLevel = 0;
    MultiGridManagement2D management(Lx, Ly, numLevel, overlapWidth, behaviorLevel);

    // Padded beyond the block's own footprint so the coarse/fine
    // refinement interface sits in guaranteed-open freestream, never
    // touching a BounceBack cell.
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
            management, new MRTdynamics<T, DESCRIPTOR>((T)1.0), behaviorLevel);

    pcout << "Lattice built with " << lattice.getNumLevels() << " levels\n";

    for (plint iLevel = 0; iLevel < lattice.getNumLevels(); ++iLevel) {
        stampPoreOnLevel(lattice, iLevel, levelDomains, blockX0, blockY0);
    }

    Array<T, 2> u0(0.02, 0.0);
    for (plint iLevel = 0; iLevel < lattice.getNumLevels(); ++iLevel) {
        MultiBlockLattice2D<T, DESCRIPTOR>& comp = lattice.getComponent(iLevel);
        initializeAtEquilibrium(comp, comp.getBoundingBox(), (T)1.0, u0);
    }
    lattice.initialize();
    for (int i = 0; i < 5; ++i) lattice.collideAndStream();

    Box2D refineBoxForWrite(blockX0, blockX0 + blockNx - 1, blockY0, blockY0 + blockNy - 1);
    for (plint iLevel = 0; iLevel < lattice.getNumLevels(); ++iLevel) {
        MultiBlockLattice2D<T, DESCRIPTOR>& comp = lattice.getComponent(iLevel);
        plint scale = 1 << iLevel;

        Box2D writeBox = comp.getBoundingBox();
        if (iLevel > 0) {
            writeBox = Box2D(refineBoxForWrite.x0 * scale, refineBoxForWrite.x1 * scale,
                              refineBoxForWrite.y0 * scale, refineBoxForWrite.y1 * scale);
        }

        VtkImageOutput2D<T> vtkOut(
            createFileName("multigrid_level" + std::to_string(iLevel), 0, 6),
            1.0 / static_cast<double>(scale));
        vtkOut.writeData<float>(*computeDensity(comp, writeBox), "density", 1.0);
    }
    pcout << "Wrote one VTK file per level to visualise the refinement structure\n";

    return 0;
}