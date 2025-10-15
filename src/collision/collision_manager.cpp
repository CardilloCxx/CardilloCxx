#include "collision_manager.hpp"

namespace cardillo::collision {

std::vector<Contact> CollisionManager::detectAll(const PhysicsSystem& sys) const {
    auto data = m_broad.collect(sys);
    auto proxies = m_broad.buildProxies(data);
    auto pairs = m_broad.makePairs(proxies);

    std::vector<Contact> contacts;
    contacts.reserve(proxies.size());

    // Dispatch pairs by collider type to the appropriate narrow-phase function
    // We batch by type combos to reuse existing narrow-phase methods.
    std::vector<std::pair<int,int>> ss; ss.reserve(pairs.size());
    std::vector<std::pair<int,int>> sp; sp.reserve(pairs.size());
    std::vector<std::pair<int,int>> so; so.reserve(pairs.size());

    for (const auto& pr : pairs) {
        const auto& A = proxies[pr.i];
        const auto& B = proxies[pr.j];
        auto tA = A.type, tB = B.type;
        // sphere-sphere
        if (tA == ColliderType::Sphere && tB == ColliderType::Sphere) {
            ss.emplace_back(A.index, B.index);
        } else if (tA == ColliderType::Sphere && tB == ColliderType::Plane) {
            sp.emplace_back(A.index, B.index);
        } else if (tA == ColliderType::Plane && tB == ColliderType::Sphere) {
            sp.emplace_back(B.index, A.index);
        } else if (tA == ColliderType::Sphere && tB == ColliderType::Obb) {
            so.emplace_back(A.index, B.index);
        } else if (tA == ColliderType::Obb && tB == ColliderType::Sphere) {
            so.emplace_back(B.index, A.index);
        }
        // other type pairs currently unsupported (plane-plane, obb-obb, plane-obb)
    }

    m_narrow.sphereSphere(data.spheres, ss, contacts);
    m_narrow.spherePlane(data.spheres, data.planes, sp, contacts);
    m_narrow.sphereObb(data.spheres, data.obbs, so, contacts);

    return contacts;
}

}
