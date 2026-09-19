// beam_hull_shape.cpp
#include "beam_hull_shape.hpp"

#include <limits>

namespace cardillo::collision {

inline std::vector<coal::Vec3s> ringFromEntityPose(const entt::registry& reg, entt::entity e, const std::vector<Vector2r>& polygon) {
    const Vector3r origin = reg.any_of<C_Position3>(e) ? reg.get<C_Position3>(e).value : Vector3r::Zero();
    const Quaternion4r q = reg.any_of<C_Orientation>(e) ? reg.get<C_Orientation>(e).value : Quaternion4r::Identity();
    std::vector<coal::Vec3s> ring;
    ring.reserve(polygon.size());
    for (const Vector2r& p2 : polygon) {
        ring.push_back(coal::Vec3s(origin + q * Vector3r((real_t)0, p2.x(), p2.y())));
    }
    return ring;
}

BeamHullShape::BeamHullShape(const C_Collider_BeamHull& collider, const entt::registry& reg){
    updateRings(collider, reg);
 }

void BeamHullShape::updateRings(const C_Collider_BeamHull& collider, const entt::registry& reg) {
    if (collider.radius > 0.0f) {
        circular_ = true;
        discA_.c = reg.get<C_Position3>(collider.endA).value;
        discB_.c = reg.get<C_Position3>(collider.endB).value;
        discA_.n = reg.get<C_Orientation>(collider.endA).rotation.col(0);
        discB_.n = reg.get<C_Orientation>(collider.endB).rotation.col(0);
        discA_.r = collider.radius;
        discB_.r = collider.radius;
    } else {
        ringA_ = ringFromEntityPose(reg, collider.endA, collider.polygon);
        ringB_ = ringFromEntityPose(reg, collider.endB, collider.polygon);
    }
    recomputeLocalAABB_();
}

void BeamHullShape::computeShapeSupport(const coal::Vec3s& dir, coal::Vec3s& support, int& hint, coal::details::ShapeSupportData& /*data*/) const {
    coal::CoalScalar bestDot = -std::numeric_limits<coal::CoalScalar>::infinity();
    support = coal::Vec3s::Zero();

    if (circular_) {
        auto discSupport = [&](const Disc& D, int discTag) {
            const coal::Vec3s perp = dir - dir.dot(D.n) * D.n;
            const coal::CoalScalar len = perp.norm();
            coal::Vec3s p;
            if (len > std::numeric_limits<coal::CoalScalar>::epsilon() * dir.norm()) {
                p = D.c + (D.r / len) * perp;
            } else {
                // dir parallel to the disc normal: every rim point is valid, so pick any perpendicular to n
                const coal::Vec3s axis = (std::abs(D.n.x()) < 0.9) ? coal::Vec3s::UnitX() : coal::Vec3s::UnitY();
                p = D.c + D.r * D.n.cross(axis).normalized();
            }
            const coal::CoalScalar s = dir.dot(D.c) + D.r * len;
            if (s > bestDot) {
                bestDot = s;
                support = p;
            }
        };

        discSupport(discA_, 0);
        discSupport(discB_, 1);
    } else {
        auto scanRing = [&](const std::vector<coal::Vec3s>& ring, int ringTag) {
            for (std::size_t i = 0; i < ring.size(); ++i) {
                const coal::CoalScalar d = dir.dot(ring[i]);
                if (d > bestDot) {
                    bestDot = d;
                    support = ring[i];
                }
            }
        };
        scanRing(ringA_, 0);
        scanRing(ringB_, 1);
    }
    hint = 0;
}

void BeamHullShape::computeLocalAABB() { recomputeLocalAABB_(); }

coal::CollisionGeometry* BeamHullShape::clone() const { return new BeamHullShape(*this); }

bool BeamHullShape::isEqual(const coal::CollisionGeometry& other) const {
    const auto* o = dynamic_cast<const BeamHullShape*>(&other);
    return o != nullptr && o->ringA_ == ringA_ && o->ringB_ == ringB_;
}

void BeamHullShape::recomputeLocalAABB_() {

    if (circular_) {
        auto half = [](const Disc& D) -> coal::Vec3s {
            return D.r * (coal::Vec3s::Ones() - D.n.cwiseProduct(D.n)).cwiseMax(0).cwiseSqrt();
        };
        const coal::Vec3s hA = half(discA_), hB = half(discB_);
        const coal::Vec3s lo = (discA_.c - hA).cwiseMin(discB_.c - hB);
        const coal::Vec3s hi = (discA_.c + hA).cwiseMax(discB_.c + hB);
        aabb_local = coal::AABB(lo, hi);
        aabb_center = aabb_local.center();
        aabb_radius = (hi - lo).norm() / 2;
        return;
    }

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
