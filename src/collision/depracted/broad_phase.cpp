#include "broad_phase.hpp"
#include <algorithm>

namespace cardillo::collision {

// Collect typed colliders
BroadPhaseData BroadPhase::collect(const PhysicsSystem& sys) const {
    BroadPhaseData out;
    const auto& reg = sys.ecs();

    // Spheres
    {
        auto view = reg.view<PhysicsSystem::C_Collidable, PhysicsSystem::C_PointMassTag, PhysicsSystem::C_Position3, PhysicsSystem::C_Radius>();
        out.spheres.reserve(view.size_hint());
        for (auto [e, pos, rad] : view.each()) out.spheres.push_back(SphereCollider{e, pos.value, rad.r});
    }

    // Planes
    {
        auto view = reg.view<PhysicsSystem::C_Collidable, PhysicsSystem::C_PlaneVisualTag, PhysicsSystem::C_Position3, PhysicsSystem::C_Plane>();
        out.planes.reserve(view.size_hint());
        for (auto [e, pos, pl] : view.each()) {
            Vector3r n = pl.normal; n.normalize();
            Vector3r upHint = pl.up;
            Vector3r u = upHint - n * (upHint.dot(n));
            if (u.norm() < (real_t)1e-8) {
                Vector3r tangent = (std::abs(n.x()) > (real_t)0.9) ? Vector3r(0,1,0) : Vector3r(1,0,0);
                u = (tangent - n * (tangent.dot(n))).normalized();
            } else u.normalize();
            Vector3r v = n.cross(u);
            out.planes.push_back(PlaneCollider{e, pos.value, n, u, v, pl.sizeX, pl.sizeY});
        }
    }

    // OBBs (cubes)
    {
        auto view = reg.view<PhysicsSystem::C_Collidable,
                             PhysicsSystem::C_CubeVisualTag,
                             PhysicsSystem::C_Position3,
                             PhysicsSystem::C_Cube,
                             PhysicsSystem::C_Orientation>();
        out.obbs.reserve(view.size_hint());
        for (auto [e, pos, cb, ori] : view.each()) {
            Matrix33r R = ori.q.toRotationMatrix();
            out.obbs.push_back(ObbCollider{e, pos.value, cb.halfExtents, R});
        }
    }
    return out;
}

// Templated aabb converters
template <typename T>
AABB BroadPhase::aabbOf(const T& t) { return aabbOf(t); }

AABB BroadPhase::aabbOf(const SphereCollider& s) {
    Vector3r r(s.radius, s.radius, s.radius);
    return {s.center - r, s.center + r};
}

AABB BroadPhase::aabbOf(const PlaneCollider& p) {
    // Finite oriented plane rectangle; computing conservative AABB from corners
    Vector3r c0 = p.point + (-p.hx)*p.u + (-p.hy)*p.v;
    Vector3r c1 = p.point + ( p.hx)*p.u + (-p.hy)*p.v;
    Vector3r c2 = p.point + ( p.hx)*p.u + ( p.hy)*p.v;
    Vector3r c3 = p.point + (-p.hx)*p.u + ( p.hy)*p.v;
    Vector3r mn = c0.cwiseMin(c1).cwiseMin(c2).cwiseMin(c3);
    Vector3r mx = c0.cwiseMax(c1).cwiseMax(c2).cwiseMax(c3);
    return {mn, mx};
}

AABB BroadPhase::aabbOf(const ObbCollider& b) {
    // Conservative AABB of OBB: center +/- |R| * he (abs rotation)
    Matrix33r absR = b.R.cwiseAbs();
    Vector3r ext = absR * b.halfExtents;
    return {b.center - ext, b.center + ext};
}

std::vector<AabbProxy> BroadPhase::buildProxies(const BroadPhaseData& data) const {
    std::vector<AabbProxy> proxies;
    proxies.reserve(data.spheres.size() + data.planes.size() + data.obbs.size());

    for (int i = 0; i < (int)data.spheres.size(); ++i)
        proxies.push_back({aabbOf(data.spheres[i]), ColliderType::Sphere, i});
    for (int i = 0; i < (int)data.planes.size(); ++i)
        proxies.push_back({aabbOf(data.planes[i]), ColliderType::Plane, i});
    for (int i = 0; i < (int)data.obbs.size(); ++i)
        proxies.push_back({aabbOf(data.obbs[i]), ColliderType::Obb, i});

    return proxies;
}

// Stable insertion sort to reuse previous order (axis.endpoints pre-sized correctly)
void BroadPhase::updateAndSortAxis(AxisData& axis, const std::vector<AabbProxy>& proxies, int axisIdx) const {
    const int N = static_cast<int>(proxies.size());
    auto& a = axis.endpoints;

    // If size changed, rebuild endpoints and do a full std::sort as a fresh baseline
    if ((int)a.size() != 2 * N) {
        a.clear();
        a.resize(2 * N);
        for (int i = 0; i < N; ++i) {
            const auto& box = proxies[i].box;
            const real_t minv = (axisIdx == 0 ? box.min.x() : (axisIdx == 1 ? box.min.y() : box.min.z()));
            const real_t maxv = (axisIdx == 0 ? box.max.x() : (axisIdx == 1 ? box.max.y() : box.max.z()));
            a[2*i + 0] = Endpoint{minv, i, true};
            a[2*i + 1] = Endpoint{maxv, i, false};
        }
        std::sort(a.begin(), a.end(), [](const Endpoint& L, const Endpoint& R){
            if (L.value == R.value) return (L.start && !R.start);
            return L.value < R.value;
        });
        return;
    }

    // Otherwise, update endpoint values in place (preserving current order)
    for (auto& ep : a) {
        const auto& box = proxies[ep.proxyIndex].box;
        ep.value = ep.start
            ? (axisIdx == 0 ? box.min.x() : (axisIdx == 1 ? box.min.y() : box.min.z()))
            : (axisIdx == 0 ? box.max.x() : (axisIdx == 1 ? box.max.y() : box.max.z()));
    }
    // Stable insertion sort for temporal coherence
    for (int i = 1; i < (int)a.size(); ++i) {
        Endpoint key = a[i];
        int j = i - 1;
        while (j >= 0) {
            const bool tie = (a[j].value == key.value) && (key.start && !a[j].start);
            if (a[j].value > key.value || tie) a[j+1] = a[j], --j; else break;
        }
        a[j+1] = key;
    }
}

void BroadPhase::sweepAxis(const AxisData& axis, std::vector<Pair>& out) const {
    out.clear();
    m_active.clear(); m_active.reserve(axis.endpoints.size());
    // Position map for O(1) removal: proxyIndex -> position in 'm_active' (or -1 if not active)
    const int N = static_cast<int>(axis.endpoints.size() / 2);
    m_pos.assign(static_cast<size_t>(N), -1);
    for (const auto& ep : axis.endpoints) {
        if (ep.start) {
            for (int j : m_active) out.push_back({ep.proxyIndex, j});
            m_pos[ep.proxyIndex] = static_cast<int>(m_active.size());
            m_active.push_back(ep.proxyIndex);
        } else {
            int p = m_pos[ep.proxyIndex];
            if (p >= 0) {
                int lastIdx = m_active.back();
                std::swap(m_active[p], m_active.back());
                m_active.pop_back();
                m_pos[lastIdx] = p;
                m_pos[ep.proxyIndex] = -1;
            }
        }
    }
}

std::vector<Pair> BroadPhase::makePairs(const std::vector<AabbProxy>& proxies) const {
    // Update and sort endpoints only for X (Y/Z used for numeric overlap checks only)
    updateAndSortAxis(m_axisX, proxies, 0);

    // Primary sweep on X (reuse buffer, no temporary)
    sweepAxis(m_axisX, m_candidates);

    // Filter with Y and Z overlaps (reuse buffer)
    m_filtered.clear(); m_filtered.reserve(m_candidates.size());
    for (const auto& pr : m_candidates) {
        const auto& Ab = proxies[pr.i].box;
        const auto& Bb = proxies[pr.j].box;
        const real_t Ay0 = Ab.min.y(), Ay1 = Ab.max.y();
        const real_t By0 = Bb.min.y(), By1 = Bb.max.y();
        if (!overlap1D(Ay0, Ay1, By0, By1)) continue;
        const real_t Az0 = Ab.min.z(), Az1 = Ab.max.z();
        const real_t Bz0 = Bb.min.z(), Bz1 = Bb.max.z();
        if (overlap1D(Az0, Az1, Bz0, Bz1)) m_filtered.push_back(pr);
    }
    return m_filtered;
}

}
