// beam_hull_shape.cpp
#include "beam_hull_shape.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include <coal/narrowphase/support_functions.h>  // details::computeSupportSetConvexHull
#include "tangent_frame.hpp"

namespace cardillo::collision {

inline std::vector<Vector3r> ringFromEntityPose(const entt::registry& reg, entt::entity e, const std::vector<Vector2r>& polygon) {
    const Vector3r origin = reg.any_of<C_Position3>(e) ? reg.get<C_Position3>(e).value : Vector3r::Zero();
    const Quaternion4r q = reg.any_of<C_Orientation>(e) ? reg.get<C_Orientation>(e).value : Quaternion4r::Identity();
    std::vector<Vector3r> ring;
    ring.reserve(polygon.size());
    for (const Vector2r& p2 : polygon) {
        ring.push_back(origin + q * Vector3r((real_t)0, p2.x(), p2.y()));
    }
    return ring;
}

BeamHullShape::BeamHullShape(const C_Collider_BeamHull& collider, const entt::registry& reg){
    updateRings(collider, reg);
 }

void BeamHullShape::updateRings(const C_Collider_BeamHull& collider, const entt::registry& reg) {
    if (collider.radiusA > 0.0f && collider.radiusB > 0.0f) {
        circular_ = true;
        ringA_.clear();
        ringB_.clear();
        discA_.c = reg.get<C_Position3>(collider.endA).value;
        discB_.c = reg.get<C_Position3>(collider.endB).value;
        discA_.n = reg.get<C_Orientation>(collider.endA).rotation.col(0);
        discB_.n = reg.get<C_Orientation>(collider.endB).rotation.col(0);
        discA_.r = collider.radiusA;
        discB_.r = collider.radiusB;
    } else {
        circular_ = false;
        discA_ = Disc{};
        discB_ = Disc{};
        ringA_ = ringFromEntityPose(reg, collider.endA, collider.polygonA);
        ringB_ = ringFromEntityPose(reg, collider.endB, collider.polygonB);
    }
    recomputeLocalAABB_();
}

Vector3r BeamHullShape::endpointSupport(const Vector3r& dir, int endpoint) const {
    if (circular_) {
        const Disc& disc = endpoint == 0 ? discA_ : discB_;
        const Vector3r perp = dir - dir.dot(disc.n) * disc.n;
        const real_t len = perp.norm();
        if (len > std::numeric_limits<real_t>::epsilon() * dir.norm()) {
            return disc.c + (disc.r / len) * perp;
        }
        const Vector3r axis = (std::abs(disc.n.x()) < 0.9) ? Vector3r::UnitX() : Vector3r::UnitY();
        return disc.c + disc.r * disc.n.cross(axis).normalized();
    }

    const auto& ring = endpoint == 0 ? ringA_ : ringB_;
    Vector3r support = Vector3r::Zero();
    real_t bestDot = -std::numeric_limits<real_t>::infinity();
    for (const auto& point : ring) {
        const real_t dot = dir.dot(point);
        if (dot > bestDot) {
            bestDot = dot;
            support = point;
        }
    }
    return support;
}

void BeamHullShape::computeShapeSupport(const Vector3r& dir, Vector3r& support, int& hint, coal::details::ShapeSupportData& /*data*/) const {
    const Vector3r supportA = endpointSupport(dir, 0);
    const Vector3r supportB = endpointSupport(dir, 1);
    support = (dir.dot(supportA) >= dir.dot(supportB)) ? supportA : supportB;
    hint = 0;
}

void BeamHullShape::computeShapeSupportSet(coal::SupportSet& supportSet, int& /*hint*/, coal::details::ShapeSupportData& data,
                                           std::size_t numSamples, real_t tol) const {
    const Vector3r n = supportSet.getNormal();

    auto forEachCandidate = [&](auto&& f) {
        if (circular_) {
            for (int e = 0; e < 2; ++e) {
                const Disc& d = e == 0 ? discA_ : discB_;

                const real_t s = std::sqrt(std::max<real_t>(0, 1 - n.dot(d.n) * n.dot(d.n)));
                if (2 * d.r * s <= tol) {
                    Vector3r u, v;
                    tangentFrameFromNormal(d.n, u, v);
                    for (std::size_t k = 0; k < numSamples; ++k) {
                        const real_t a = 2 * M_PI * real_t(k) / real_t(numSamples);
                        f(d.c + d.r * (std::cos(a) * u + std::sin(a) * v));
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

    real_t best = -std::numeric_limits<real_t>::infinity();
    forEachCandidate([&](const Vector3r& p) { best = std::max(best, n.dot(p)); });

    auto& poly = data.polygon;
    poly.clear();
    forEachCandidate([&](const Vector3r& p) {
        if (best - n.dot(p) <= tol) poly.push_back(supportSet.tf.inverseTransform(p).head<2>());
    });
    coal::details::computeSupportSetConvexHull(poly, supportSet.points());
}

void BeamHullShape::computeLocalAABB() { recomputeLocalAABB_(); }

coal::CollisionGeometry* BeamHullShape::clone() const { return new BeamHullShape(*this); }

bool BeamHullShape::isEqual(const coal::CollisionGeometry& other) const {
    const auto* o = dynamic_cast<const BeamHullShape*>(&other);
    if (o == nullptr) return false;

    COAL_EQUAL_OPERATOR_CHECK(ShapeBase::isEqual(other));
    COAL_EQUAL_OPERATOR_CHECK(circular_ == o->circular_);
    if (circular_) {
        COAL_EQUAL_OPERATOR_CHECK(discA_ == o->discA_);
        COAL_EQUAL_OPERATOR_CHECK(discB_ == o->discB_);
    } else {
        COAL_EQUAL_OPERATOR_CHECK(ringA_ == o->ringA_);
        COAL_EQUAL_OPERATOR_CHECK(ringB_ == o->ringB_);
    }
    return true;
}

void BeamHullShape::recomputeLocalAABB_() {

    if (circular_) {
        auto half = [](const Disc& D) -> Vector3r {
            return D.r * (Vector3r::Ones() - D.n.cwiseProduct(D.n)).cwiseMax(0).cwiseSqrt();
        };
        const Vector3r hA = half(discA_), hB = half(discB_);
        const Vector3r lo = (discA_.c - hA).cwiseMin(discB_.c - hB);
        const Vector3r hi = (discA_.c + hA).cwiseMax(discB_.c + hB);
        aabb_local = coal::AABB(lo, hi);
        aabb_center = aabb_local.center();
        aabb_radius = (hi - lo).norm() / 2;
        return;
    }

    if (ringA_.empty() && ringB_.empty()) {
        aabb_local = coal::AABB(Vector3r::Zero());
        aabb_center = Vector3r::Zero();
        aabb_radius = 0;
        return;
    }
    Vector3r lo = ringA_.empty() ? ringB_.front() : ringA_.front();
    Vector3r hi = lo;
    auto absorb = [&](const std::vector<Vector3r>& ring) {
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