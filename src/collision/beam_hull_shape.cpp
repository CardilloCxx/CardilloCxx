// beam_hull_shape.cpp
#include "beam_hull_shape.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include <coal/narrowphase/support_functions.h>  // details::computeSupportSetConvexHull

namespace cardillo::collision {

namespace {
// Orthonormal basis (u, v) of the plane orthogonal to the unit vector n.
inline void orthonormalBasis(const coal::Vec3s& n, coal::Vec3s& u, coal::Vec3s& v) {
    const coal::Vec3s axis = (std::abs(n.x()) < 0.9) ? coal::Vec3s::UnitX() : coal::Vec3s::UnitY();
    u = n.cross(axis).normalized();
    v = n.cross(u);
}
}  // namespace

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

coal::Vec3s BeamHullShape::endpointSupport(const coal::Vec3s& dir, int endpoint) const {
    if (circular_) {
        const Disc& disc = endpoint == 0 ? discA_ : discB_;
        const coal::Vec3s perp = dir - dir.dot(disc.n) * disc.n;
        const coal::CoalScalar len = perp.norm();
        if (len > std::numeric_limits<coal::CoalScalar>::epsilon() * dir.norm()) {
            return disc.c + (disc.r / len) * perp;
        }
        const coal::Vec3s axis = (std::abs(disc.n.x()) < 0.9) ? coal::Vec3s::UnitX() : coal::Vec3s::UnitY();
        return disc.c + disc.r * disc.n.cross(axis).normalized();
    }

    const auto& ring = endpoint == 0 ? ringA_ : ringB_;
    coal::Vec3s support = coal::Vec3s::Zero();
    coal::CoalScalar bestDot = -std::numeric_limits<coal::CoalScalar>::infinity();
    for (const auto& point : ring) {
        const coal::CoalScalar dot = dir.dot(point);
        if (dot > bestDot) {
            bestDot = dot;
            support = point;
        }
    }
    return support;
}

void BeamHullShape::computeShapeSupport(const coal::Vec3s& dir, coal::Vec3s& support, int& hint, coal::details::ShapeSupportData& /*data*/) const {
    const coal::Vec3s supportA = endpointSupport(dir, 0);
    const coal::Vec3s supportB = endpointSupport(dir, 1);
    support = (dir.dot(supportA) >= dir.dot(supportB)) ? supportA : supportB;
    hint = 0;
}

void BeamHullShape::computeShapeSupportSet(coal::SupportSet& supportSet, int& /*hint*/, coal::details::ShapeSupportData& data,
                                           std::size_t numSamples, coal::CoalScalar tol) const {
    // Support direction in this shape's local frame; already flipped when the hull is the 2nd shape.
    const coal::Vec3s n = supportSet.getNormal();

    // Visits every candidate point of the support set. The slab test needs two passes (extreme value
    // first, then collect), so candidates are generated on the fly instead of being stored.
    auto forEachCandidate = [&](auto&& f) {
        if (circular_) {
            for (int e = 0; e < 2; ++e) {
                const Disc& d = e == 0 ? discA_ : discB_;
                // Extent of the rim along n is 2*r*s. If the whole rim lies within tol of the extreme
                // plane, the disc itself is the contact face (end cap): return a fixed-frame polygon of
                // rim samples so the patch does not rotate from frame to frame. Otherwise only the exact
                // extreme rim point can touch (beam lying on its side): together, the two discs give
                // the contact segment.
                const coal::CoalScalar s = std::sqrt(std::max<coal::CoalScalar>(0, 1 - n.dot(d.n) * n.dot(d.n)));
                if (2 * d.r * s <= tol) {
                    coal::Vec3s u, v;
                    orthonormalBasis(d.n, u, v);
                    for (std::size_t k = 0; k < numSamples; ++k) {
                        const coal::CoalScalar a = 2 * M_PI * coal::CoalScalar(k) / coal::CoalScalar(numSamples);
                        f(coal::Vec3s(d.c + d.r * (std::cos(a) * u + std::sin(a) * v)));
                    }
                } else {
                    f(endpointSupport(n, e));
                }
            }
        } else {
            for (const auto& p : ringA_) f(p);
            for (const auto& p : ringB_) f(p);
        }
    };

    coal::CoalScalar best = -std::numeric_limits<coal::CoalScalar>::infinity();
    forEachCandidate([&](const coal::Vec3s& p) { best = std::max(best, n.dot(p)); });

    // Keep the points within tol of the extreme plane, as Coal does for Box and ConvexBase, and
    // hand their 2D convex hull (in the patch frame) to the contact patch solver.
    auto& poly = data.polygon;
    poly.clear();
    forEachCandidate([&](const coal::Vec3s& p) {
        if (best - n.dot(p) <= tol) poly.push_back(supportSet.tf.inverseTransform(p).head<2>());
    });
    coal::details::computeSupportSetConvexHull(poly, supportSet.points());
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