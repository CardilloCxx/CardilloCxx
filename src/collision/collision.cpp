#include "collision.hpp"
#include <algorithm>
#include <cmath>

namespace cardillo::collision {

using PS = cardillo::PhysicsSystem;

// Generic fallback: unsupported pair -> no contact
template <typename A, typename B>
static std::optional<Contact> detect(const A&, const B&) {
    return std::nullopt;
}

// Static dispatch detect() overloads
static std::optional<Contact> detect(const SphereCollider& a, const SphereCollider& b) {
    Vector3r d = b.center - a.center;
    real_t dist = d.norm();
    real_t rsum = a.radius + b.radius;
    if (dist < rsum) {
        Vector3r n = (dist > (real_t)1e-12) ? (d / dist) : Vector3r(1,0,0);
        real_t penetration = rsum - dist;
        Vector3r contact = a.center + n * (a.radius - (real_t)0.5 * penetration);
        return Contact{a.e, b.e, contact, n, penetration};
    }
    return std::nullopt;
}

static std::optional<Contact> detect(const PlaneCollider& p, const SphereCollider& s) {
    // distance center sphere and origin plane
    Vector3r rel = s.center - p.point;

    // distance along normal
    const Vector3r n = p.normal;
    real_t d = rel.dot(n);

    if (d < -s.radius) return std::nullopt; // sphere behind plane

    // coordinates along plane axes
    real_t u = rel.dot(p.u);
    real_t v = rel.dot(p.v);
    real_t cu = std::max(-p.hx, std::min(u, p.hx));
    real_t cv = std::max(-p.hy, std::min(v, p.hy));

    if ((u < -p.hx) || (p.hx < u) || (v < -p.hy) || (p.hy < v))
        return std::nullopt;

    real_t penetration = d - s.radius;
    if (penetration > 0) return std::nullopt;
        
    Vector3r contact = p.point + cu * p.u + cv * p.v; // closest point on finite plane
    return Contact{s.e, p.e, contact, p.normal, penetration};
}

static std::optional<Contact> detect(const SphereCollider& s, const PlaneCollider& p) {
    if (auto c = detect(p, s)) {
        // swap so normal points from p->s consistently
        Contact out = *c;
        out.a = p.e; out.b = s.e; out.normal = out.normal; // swap normal?
        return out;
    }
    return std::nullopt;
}

// static std::optional<Contact> detect(const SphereCollider& s, const PlaneCollider& p) {
//     // distance along normal
//     const Vector3r n = p.normal;
//     real_t d = (s.center - p.point).dot(n);
//     if (std::abs(d) > s.radius) return std::nullopt; // no intersection with plane slab (infinite plane)

//     // project onto plane and clamp within rectangle bounds
//     Vector3r rel = s.center - p.point - d * n; // in-plane vector
//     // coordinates along plane axes
//     real_t u = rel.dot(p.u);
//     real_t v = rel.dot(p.v);
//     real_t cu = std::max(-p.hx, std::min(p.hx, u));
//     real_t cv = std::max(-p.hy, std::min(p.hy, v));
//     Vector3r closest = p.point + cu * p.u + cv * p.v; // closest point on finite plane
//     Vector3r diff = s.center - closest; // vector from plane to sphere
//     real_t dist = diff.norm();
//     if (dist <= s.radius) {
//         real_t penetration = s.radius - dist;
//         Vector3r nout;
//         if (dist > (real_t)1e-12) {
//             // normal should point sphere -> plane
//             nout = -diff / dist;
//         } else {
//             // use sign of signed distance d to choose direction consistently (sphere->plane)
//             nout = (d >= 0) ? -n : n;
//         }
//         Vector3r contact = closest; // approximate contact point at plane
//         return Contact{s.e, p.e, contact, nout, penetration};
//     }
//     return std::nullopt;
// }

// static std::optional<Contact> detect(const PlaneCollider& p, const SphereCollider& s) {
//     if (auto c = detect(s, p)) {
//         // swap so normal points from p->s consistently
//         Contact out = *c;
//         out.a = p.e; out.b = s.e; out.normal = -out.normal;
//         return out;
//     }
//     return std::nullopt;
// }

static std::optional<Contact> detect(const SphereCollider& s, const AabbCollider& b) {
    Vector3r bmin = b.center - b.halfExtents;
    Vector3r bmax = b.center + b.halfExtents;
    Vector3r closest = s.center.cwiseMax(bmin).cwiseMin(bmax);
    Vector3r d = s.center - closest; // from box to sphere
    real_t dist2 = d.squaredNorm();
    if (dist2 <= s.radius * s.radius) {
        real_t dist = std::sqrt((double)dist2);
        Vector3r n;
        if (dist > (real_t)1e-12) {
            // normal sphere->box is opposite to d
            n = -d / dist;
        } else {
            // choose axis toward nearest face (sphere->box)
            Vector3r delta = s.center - b.center;
            Vector3r faceDist(
                b.halfExtents.x() - std::abs(delta.x()),
                b.halfExtents.y() - std::abs(delta.y()),
                b.halfExtents.z() - std::abs(delta.z())
            );
            int axis = 0;
            if (faceDist.y() < faceDist[axis]) axis = 1;
            if (faceDist.z() < faceDist[axis]) axis = 2;
            n = Vector3r::Zero();
            real_t sgn = (delta[axis] >= 0) ? 1 : -1;
            n[axis] = -sgn; // sphere->box direction is opposite to outward face normal
        }
        real_t penetration = s.radius - dist;
        Vector3r contact = closest;
        return Contact{s.e, b.e, contact, n, penetration};
    }
    return std::nullopt;
}

static std::optional<Contact> detect(const AabbCollider& b, const SphereCollider& s) {
    if (auto c = detect(s, b)) {
        Contact out = *c; out.a = b.e; out.b = s.e; out.normal = -out.normal; return out;
    }
    return std::nullopt;
}

// OBB helper: compute closest point on OBB to a point
static Vector3r closestPointOnObb(const Vector3r& p, const Vector3r& c, const Matrix33r& R, const Vector3r& he) {
    Vector3r d = p - c;
    Vector3r q = c;
    for (int i = 0; i < 3; ++i) {
        Vector3r axis = R.col(i);
        real_t dist = d.dot(axis);
        real_t clamped = std::max(-he[i], std::min(he[i], dist));
        q += clamped * axis;
    }
    return q;
}

static std::optional<Contact> detect(const SphereCollider& s, const ObbCollider& b) {
    Vector3r closest = closestPointOnObb(s.center, b.center, b.R, b.halfExtents);
    Vector3r d = s.center - closest; // from box to sphere
    real_t dist2 = d.squaredNorm();
    if (dist2 <= s.radius * s.radius) {
        real_t dist = std::sqrt((double)dist2);
        Vector3r n;
        if (dist > (real_t)1e-12) {
            // normal sphere->box is opposite to d
            n = -d / dist;
        } else {
            // choose nearest face normal in OBB local space, then map to world; keep sphere->box direction
            Vector3r local = b.R.transpose() * (s.center - b.center);
            Vector3r faceDist(
                b.halfExtents.x() - std::abs(local.x()),
                b.halfExtents.y() - std::abs(local.y()),
                b.halfExtents.z() - std::abs(local.z())
            );
            int axis = 0;
            if (faceDist.y() < faceDist[axis]) axis = 1;
            if (faceDist.z() < faceDist[axis]) axis = 2;
            real_t sgn = (local[axis] >= 0) ? 1 : -1;
            n = -sgn * b.R.col(axis); // opposite of outward normal for sphere->box direction
        }
        real_t penetration = s.radius - dist;
        Vector3r contact = closest;
        return Contact{s.e, b.e, contact, n, penetration}; 
    }
    return std::nullopt;
}

static std::optional<Contact> detect(const ObbCollider& b, const SphereCollider& s) {
    if (auto c = detect(s, b)) { Contact out = *c; out.a = b.e; out.b = s.e; out.normal = -out.normal; return out; }
    return std::nullopt;
}

// Robust plane basis construction: try using upHint, fallback to a stable axis
static inline void makePlaneBasis(const Vector3r& n, const Vector3r& upHint, Vector3r& outU, Vector3r& outV) {
    Vector3r upProj = upHint - n * (upHint.dot(n));
    if (upProj.norm() > (real_t)1e-8) {
        outU = upProj.normalized();
    } else {
        Vector3r tangent = (std::abs(n.x()) > (real_t)0.9) ? Vector3r(0,1,0) : Vector3r(1,0,0);
        outU = (tangent - n * (tangent.dot(n))).normalized();
    }
    outV = n.cross(outU);
}

std::vector<Contact> detectAll(const PhysicsSystem& sys) {
    std::vector<Contact> contacts;
    const auto& reg = sys.ecs();

    // Collect generic colliders
    std::vector<Collider> colliders;
    // Spheres (point masses)
    {
        auto view = reg.view<PS::C_Collidable, PS::C_PointMassTag, PS::C_Position3, PS::C_Radius>();
        for (auto [e, pos, rad] : view.each()) {
            colliders.emplace_back(SphereCollider{e, pos.value, rad.r});
        }
    }
    // Planes (finite, oriented)
    {
        auto view = reg.view<PS::C_Collidable, PS::C_PlaneVisualTag, PS::C_Position3, PS::C_Plane>();
        for (auto [e, pos, pl] : view.each()) {
            Vector3r n = pl.normal; n.normalize();
            Vector3r u, v;
            makePlaneBasis(n, pl.up, u, v);
            colliders.emplace_back(PlaneCollider{e, pos.value, n, u, v, pl.sizeX, pl.sizeY});
        }
    }
    // Cubes -> OBB (with rotation); also keep AABB as fallback if R ~ I
    {
        auto view = reg.view<PS::C_Collidable, PS::C_CubeVisualTag, PS::C_Position3, PS::C_Cube>();
        for (auto [e, pos, cb] : view.each()) {
            // If rotation is identity within tolerance, we could use AABB, but variant can hold OBB directly
            colliders.emplace_back(ObbCollider{e, pos.value, cb.halfExtents, cb.R});
        }
    }

    // Naive narrow phase: test all unordered pairs using variant visitation
    for (size_t i = 0; i < colliders.size(); ++i) {
        for (size_t j = i + 1; j < colliders.size(); ++j) {
            const Collider& A = colliders[i];
            const Collider& B = colliders[j];
            std::optional<Contact> c;
            std::visit([&](const auto& a){
                std::visit([&](const auto& b){ c = detect(a,b); }, B);
            }, A);
            if (c) contacts.push_back(*c);
        }
    }

    return contacts;
}

}
