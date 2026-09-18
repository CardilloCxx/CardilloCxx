// beam_hull_shape.cpp
#include "beam_hull_shape.hpp"

#include <limits>

namespace cardillo::collision {

BeamHullShape::BeamHullShape(std::vector<coal::Vec3s> ringA, std::vector<coal::Vec3s> ringB) : ringA_(std::move(ringA)), ringB_(std::move(ringB)) {
    recomputeLocalAABB_();
}

void BeamHullShape::updateRings(std::vector<coal::Vec3s> ringA, std::vector<coal::Vec3s> ringB) {
    ringA_ = std::move(ringA);
    ringB_ = std::move(ringB);
    recomputeLocalAABB_();
}

void BeamHullShape::computeShapeSupport(const coal::Vec3s& dir, coal::Vec3s& support, int& hint, coal::details::ShapeSupportData& /*data*/) const {
    coal::CoalScalar bestDot = -std::numeric_limits<coal::CoalScalar>::infinity();
    int bestHint = 0;
    support = coal::Vec3s::Zero();

    auto scanRing = [&](const std::vector<coal::Vec3s>& ring, int ringTag) {
        for (std::size_t i = 0; i < ring.size(); ++i) {
            const coal::CoalScalar d = dir.dot(ring[i]);
            if (d > bestDot) {
                bestDot = d;
                support = ring[i];
                bestHint = ringTag * 1000000 + static_cast<int>(i);
            }
        }
    };
    scanRing(ringA_, 0);
    scanRing(ringB_, 1);
    hint = bestHint;  // not load-bearing for a linear scan; kept for API parity with other shapes
}

void BeamHullShape::computeLocalAABB() { recomputeLocalAABB_(); }

coal::CollisionGeometry* BeamHullShape::clone() const { return new BeamHullShape(*this); }

bool BeamHullShape::isEqual(const coal::CollisionGeometry& other) const {
    const auto* o = dynamic_cast<const BeamHullShape*>(&other);
    return o != nullptr && o->ringA_ == ringA_ && o->ringB_ == ringB_;
}

void BeamHullShape::recomputeLocalAABB_() {
    // ringA_/ringB_ are never both empty in practice (a beam always has a real polygon per end);
    // guard anyway so a malformed C_RB_Beam degrades to a zero-size AABB instead of crashing.
    if (ringA_.empty() && ringB_.empty()) {
        aabb_local = coal::AABB(coal::Vec3s::Zero());
        aabb_center = coal::Vec3s::Zero();
        aabb_radius = 0;
        return;
    }
    coal::Vec3s lo = ringA_.empty() ? ringB_.front() : ringA_.front();
    coal::Vec3s hi = lo;
    auto absorb = [&](const std::vector<coal::Vec3s>& ring) {
        for (const auto& p : ring) {
            lo = lo.cwiseMin(p);
            hi = hi.cwiseMax(p);
        }
    };
    absorb(ringA_);
    absorb(ringB_);
    aabb_local = coal::AABB(lo, hi);
    aabb_center = aabb_local.center();
    aabb_radius = (hi - lo).norm() / 2;
}

}  // namespace cardillo::collision
