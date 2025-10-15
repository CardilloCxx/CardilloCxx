#include "narrow_phase.hpp"
#include <algorithm>
#include <cmath>

namespace cardillo::collision {

void NarrowPhase::sphereSphere(const std::vector<SphereCollider>& S,
                               const std::vector<std::pair<int,int>>& pairs,
                               std::vector<Contact>& out) const {
    for (auto [ia, ib] : pairs) {
        const auto& a = S[ia];
        const auto& b = S[ib];
        Vector3r d = b.center - a.center;
        real_t rsum = a.radius + b.radius;
        real_t rsum2 = rsum * rsum;
        real_t dist2 = d.squaredNorm();
        if (dist2 < rsum2) {
            real_t dist = std::sqrt((double)dist2);
            Vector3r n = (dist > (real_t)1e-12) ? (d / dist) : Vector3r(1,0,0);
            real_t penetration = rsum - dist;
            Vector3r contact = a.center + n * (a.radius - (real_t)0.5 * penetration);
            MatrixXXr wA(1,3); wA << -n[0], -n[1], -n[2];
            MatrixXXr wB(1,3); wB <<  n[0],  n[1],  n[2];
            out.push_back(Contact{a.e, b.e, contact, n, wA, wB, penetration});
        }
    }
}

void NarrowPhase::spherePlane(const std::vector<SphereCollider>& spheres,
                              const std::vector<PlaneCollider>& planes,
                              const std::vector<std::pair<int,int>>& pairs,
                              std::vector<Contact>& out) const {
    for (auto [is, ip] : pairs) {
        const auto& s = spheres[is];
        const auto& p = planes[ip];
        const Vector3r& n = p.normal;
        real_t d = (s.center - p.point).dot(n);
        if (std::abs(d) > s.radius) continue; // outside plane slab
        Vector3r rel = s.center - p.point - d * n;
        real_t u = rel.dot(p.u);
        real_t v = rel.dot(p.v);
        real_t cu = std::max(-p.hx, std::min(p.hx, u));
        real_t cv = std::max(-p.hy, std::min(p.hy, v));
        Vector3r closest = p.point + cu * p.u + cv * p.v;
        Vector3r diff = s.center - closest; // plane->sphere
        real_t dist = diff.norm();
        if (dist <= s.radius) {
            real_t penetration = s.radius - dist;
            Vector3r nout;
            if (dist > (real_t)1e-12) nout = -diff / dist; // sphere->plane
            else nout = (d >= 0) ? -n : n;
            MatrixXXr wA(1,3); wA << -nout[0], -nout[1], -nout[2];
            MatrixXXr wB(1,3); wB <<  nout[0],  nout[1],  nout[2];
            out.push_back(Contact{s.e, p.e, closest, nout, wA, wB, penetration});
        }
    }
}

static inline Vector3r closestPointOnObb(const Vector3r& p, const Vector3r& c, const Matrix33r& R, const Vector3r& he) {
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

void NarrowPhase::sphereObb(const std::vector<SphereCollider>& spheres,
                            const std::vector<ObbCollider>& obbs,
                            const std::vector<std::pair<int,int>>& pairs,
                            std::vector<Contact>& out) const {
    for (auto [is, ib] : pairs) {
        const auto& s = spheres[is];
        const auto& b = obbs[ib];
        Vector3r closest = closestPointOnObb(s.center, b.center, b.R, b.halfExtents);
        Vector3r d = s.center - closest; // box->sphere
        real_t dist2 = d.squaredNorm();
        if (dist2 <= s.radius * s.radius) {
            real_t dist = std::sqrt((double)dist2);
            Vector3r n;
            if (dist > (real_t)1e-12) n = -d / dist; // sphere->box
            else {
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
                n = -sgn * b.R.col(axis);
            }
            real_t penetration = s.radius - ((dist2 > 0) ? std::sqrt((double)dist2) : (real_t)0);
            MatrixXXr wA(1,3); wA << -n[0], -n[1], -n[2];
            MatrixXXr wB(1,3); wB <<  n[0],  n[1],  n[2];
            out.push_back(Contact{s.e, b.e, closest, n, wA, wB, penetration});
        }
    }
}

}
