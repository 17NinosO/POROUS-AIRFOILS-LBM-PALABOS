// PorousChannel/pore_mask_domain.hpp
#pragma once
#include "domain.hpp"

// Answers, for any lattice node (iX, iY): is this solid?
// The (row, col) <-> (iX, iY) swap lives HERE and nowhere else,
// so if it's wrong, it's wrong in exactly one place we can check.
class PoreMaskDomain {
public:
    explicit PoreMaskDomain(pore_geometry::DomainMask mask) : mask_(std::move(mask)) {}

    bool isSolid(int iX, int iY) const {
        return !mask_.at(/*row=*/iY, /*col=*/iX);
    }

    int nx() const { return mask_.nx(); }
    int ny() const { return mask_.ny(); }

private:
    pore_geometry::DomainMask mask_;
};