// COAL-based collision manager implementation (persistent, broadphase-backed)
#include "collision_coal.hpp"
#include "types.hpp"
#include "../physics/physics_system.hpp"
#include "../config/config.hpp"

#include <coal/math/transform.h>
#include <coal/shape/geometric_shapes.h>
#include <coal/collision.h>
#include <coal/collision_data.h>
#include <coal/contact_patch.h>
#include <coal/broadphase/broadphase.h>
#include <coal/broadphase/default_broadphase_callbacks.h>

#include <memory>
#include <optional>
#include <cmath>
#include <limits>
#include "tangent_frame.hpp"
#include "transform_utils.hpp"

namespace cardillo::collision {
CollisionCoal::~CollisionCoal() = default;

namespace {
inline coal::Quatf toCoalQuat(const Quaternion4r& q)
{
    return coal::Quatf(q.w(), q.x(), q.y(), q.z());
}
inline coal::Vec3s toCoalVec3(const Vector3r& v)
{
    return coal::Vec3s(v.x(), v.y(), v.z());
}
inline coal::Transform3s makeTfFromEcs(const entt::registry& reg, entt::entity e)
{
    coal::Transform3s X; X.setIdentity();
    Vector3r x = Vector3r::Zero();
    Quaternion4r q = Quaternion4r::Identity();
    if (reg.any_of<cardillo::PhysicsSystem::C_Position3>(e))
        x = reg.get<cardillo::PhysicsSystem::C_Position3>(e).value;
    if (reg.any_of<cardillo::PhysicsSystem::C_Orientation>(e))
        q = reg.get<cardillo::PhysicsSystem::C_Orientation>(e).q;
    q.normalize();
    X.setTranslation(toCoalVec3(x));
    X.setQuatRotation(toCoalQuat(q));
    return X;
}

// Helper to populate a Contact from world data and entities
static inline real_t combineFrictionMu(const entt::registry& reg, entt::entity a, entt::entity b, const std::string& method)
{
    auto getMu = [&](entt::entity e)->real_t{
        if (reg.any_of<cardillo::PhysicsSystem::C_Friction>(e)) {
            return std::max<real_t>((real_t)0, reg.get<cardillo::PhysicsSystem::C_Friction>(e).mu);
        }
        return (real_t)0;
    };
    const real_t muA = getMu(a);
    const real_t muB = getMu(b);
    if (method == "arith") return (muA + muB) * (real_t)0.5;
    if (method == "geom")  return (real_t)std::sqrt(std::max<real_t>((real_t)0, muA * muB));
    return std::min(muA, muB); // default: min
}

// Helper to populate a Contact from world data and entities
inline Contact makeContact(const entt::registry& reg,
                           entt::entity ea,
                           entt::entity eb,
                           const Vector3r& p1W,
                           const Vector3r& p2W,
                           const Vector3r& nW,
                           real_t depth,
                           const std::string& frictionCombine)
{
    // Canonicalize entity order so contact orientation is stable across frames
    entt::entity a_can = ea;
    entt::entity b_can = eb;
    Vector3r nW_can = nW;
    if (entt::to_integral(a_can) > entt::to_integral(b_can)) {
        std::swap(a_can, b_can);
        nW_can = -nW_can; // keep normal pointing from A->B
    }
    Contact c{}; c.a = a_can; c.b = b_can;
    // Precompute transforms for both bodies once
    const cardillo::collision::BodyXform XA = cardillo::collision::BodyXform::fromEcs(reg, c.a);
    const cardillo::collision::BodyXform XB = cardillo::collision::BodyXform::fromEcs(reg, c.b);
    // Normalize normal and construct stable tangent frame
    Vector3r nN = nW_can;
    real_t nlen = nN.norm();
    if (nlen > (real_t)0) nN /= nlen;
    c.normal = nN; // convention: from A -> B
    cardillo::collision::tangentFrameFromNormal(nN, c.tangent1, c.tangent2);
    c.point = (p1W + p2W) * (real_t)0.5;
    c.penetration = std::max<real_t>(0.0, depth);
    c.pointA_body  = XA.worldPointToBody(c.point);
    c.normalA_body = XA.worldVecToBody(c.normal);
    c.pointB_body  = XB.worldPointToBody(c.point);
    c.normalB_body = XB.worldVecToBody(c.normal);
    // Tangents in body frames
    c.tangent1A_body = XA.worldVecToBody(c.tangent1);
    c.tangent2A_body = XA.worldVecToBody(c.tangent2);
    c.tangent1B_body = XB.worldVecToBody(c.tangent1);
    c.tangent2B_body = XB.worldVecToBody(c.tangent2);
    // Friction coefficient combination (entity components optional)
    c.friction_mu = combineFrictionMu(reg, c.a, c.b, frictionCombine);
    return c;
}

// Insert a contact into a ContactMap bucketed by (a,b)
inline void addContactToMap(ContactMap& cmap, const Contact& c) {
    const ContactPairKey key = ContactPairKey::make(c.a, c.b);
    auto it = cmap.find(key);
    if (it == cmap.end()) {
        cmap.emplace(key, ContactList{c});
    } else {
        it->second.emplace_back(c);
    }
}

// For each current contact, find the closest previous contact (by world-point distance).
// If the closest is within maxDist, write its index to prev_pair_contact_idx; otherwise set -1.
inline void matchContactsForPair(const ContactList& prev, ContactList& curr, real_t maxDist) {
    if (curr.empty()) return;
    const real_t max2 = maxDist * maxDist;
    for (std::size_t i = 0; i < curr.size(); ++i) {
        int best = -1;
        real_t best2 = std::numeric_limits<real_t>::infinity();
        const Vector3r& p = curr[i].point;
        for (std::size_t j = 0; j < prev.size(); ++j) {
            const Vector3r d = p - prev[j].point;
            const real_t d2 = d.squaredNorm();
            if (d2 < best2) { best2 = d2; best = (int)j; }
        }
        if (best >= 0 && best2 <= max2) {
            // Record the previous contact's global index
            curr[i].prev_global_out_index = prev[best].global_out_index;
        } else {
            curr[i].prev_global_out_index = -1;
        }
    }
}

// Append contacts for a pair into a map: prefer patch expansion; fall back to raw contacts
inline void appendContactsFromPair(const entt::registry& reg,
                                   entt::entity ea,
                                   entt::entity eb,
                                   const coal::CollisionResult& cres,
                                   const coal::ContactPatchResult& patch_res,
                                   ContactMap& outMap,
                                   const std::string& frictionCombine)
{
    if (patch_res.numContactPatches() > 0) {
        for (std::size_t ip = 0; ip < patch_res.numContactPatches(); ++ip) {
            const coal::ContactPatch& patch = patch_res.getContactPatch(ip);
            const Vector3r nW(patch.getNormal().x(), patch.getNormal().y(), patch.getNormal().z());
            const real_t depth = (real_t)patch.penetration_depth;
            const std::size_t m = patch.size();
            for (std::size_t iv = 0; iv < m; ++iv) {
                const auto p1Wc = patch.getPointShape1(iv);
                const auto p2Wc = patch.getPointShape2(iv);
                const Vector3r p1W(p1Wc.x(), p1Wc.y(), p1Wc.z());
                const Vector3r p2W(p2Wc.x(), p2Wc.y(), p2Wc.z());
                addContactToMap(outMap, makeContact(reg, ea, eb, p1W, p2W, nW, depth, frictionCombine));
            }
        }
        return;
    }

    // Fallback: emit raw contacts from narrowphase result
    for (int k = 0; k < cres.numContacts(); ++k) {
        const coal::Contact c0 = cres.getContact(k);
        const Vector3r nW(c0.normal.x(), c0.normal.y(), c0.normal.z());
        const Vector3r p1W(c0.nearest_points[0].x(), c0.nearest_points[0].y(), c0.nearest_points[0].z());
        const Vector3r p2W(c0.nearest_points[1].x(), c0.nearest_points[1].y(), c0.nearest_points[1].z());
        addContactToMap(outMap, makeContact(reg, ea, eb, p1W, p2W, nW, (real_t)c0.penetration_depth, frictionCombine));
    }
}
}

void CollisionCoal::ensureBroadphaseFromConfig_() {
    if (m_broadphase) return;
    // Select based on config
    std::string bp = m_sys ? m_sys->config().collision_broadphase : std::string("dynamic_aabb");
    auto toLower = [](std::string s){ for (auto& c : s) c = (char)std::tolower((unsigned char)c); return s; };
    bp = toLower(bp);
    if (bp == "dynamic_aabb" || bp == "daabb" || bp == "default") {
        m_broadphase = std::make_unique<coal::DynamicAABBTreeCollisionManager>();
    } else if (bp == "dynamic_aabb_array" || bp == "daabb_array") {
        m_broadphase = std::make_unique<coal::DynamicAABBTreeArrayCollisionManager>();
    } else if (bp == "naive" || bp == "bruteforce") {
        m_broadphase = std::make_unique<coal::NaiveCollisionManager>();
    } else if (bp == "sap") {
        m_broadphase = std::make_unique<coal::SaPCollisionManager>();
    } else if (bp == "ssap") {
        m_broadphase = std::make_unique<coal::SSaPCollisionManager>();
    } else {
        // Fallback
        m_broadphase = std::make_unique<coal::DynamicAABBTreeCollisionManager>();
    }
}

CollisionCoal::ColliderKind CollisionCoal::inferKind_(entt::entity e) const {
    const auto& reg = m_sys->ecs();
    if (reg.any_of<PhysicsSystem::C_RB_Cube>(e)) return ColliderKind::Box;
    if (reg.any_of<PhysicsSystem::C_RB_Plane>(e)) return ColliderKind::Halfspace;
    if (reg.any_of<PhysicsSystem::C_PointMassTag>(e) && reg.any_of<PhysicsSystem::C_Radius>(e)) return ColliderKind::Sphere;
    throw std::runtime_error("CollisionCoal: unsupported collider entity; add appropriate tag/geometry.");
}

std::shared_ptr<coal::CollisionGeometry> CollisionCoal::makeGeometryFor_(ColliderKind kind, entt::entity e) const {
    const auto& reg = m_sys->ecs();
    switch (kind) {
        case ColliderKind::Box: {
            const auto& he = reg.get<PhysicsSystem::C_RB_Cube>(e).halfExtents;
            return std::make_shared<coal::Box>(he.x()*2.0, he.y()*2.0, he.z()*2.0);
        }
        case ColliderKind::Halfspace: {
            const auto& plane = reg.get<PhysicsSystem::C_RB_Plane>(e);
            Vector3r n = plane.normal.normalized();
            const auto& x = reg.get<PhysicsSystem::C_Position3>(e).value;
            real_t d = n.dot(x);
            return std::make_shared<coal::Halfspace>(toCoalVec3(n), (coal::CoalScalar)d);
        }
        case ColliderKind::Sphere: {
            const auto r = reg.get<PhysicsSystem::C_Radius>(e).r;
            return std::make_shared<coal::Sphere>((coal::CoalScalar)r);
        }
    }
    throw std::runtime_error("CollisionCoal: unknown collider kind");
}

void CollisionCoal::rebuild() {
    if (!m_sys) return;
    clear();
    ensureBroadphaseFromConfig_();
    const auto& reg = m_sys->ecs();

    // Gather all collidable entities in a single pass
    m_entities.clear();
    m_index_from_entity.clear();
    {
        auto view = reg.view<PhysicsSystem::C_Collidable, PhysicsSystem::C_Position3>();
        m_entities.reserve(view.size_hint());
        for (auto e : view) {
            m_index_from_entity.emplace(entt::to_integral(e), m_entities.size());
            m_entities.push_back(static_cast<entt::entity>(e));
        }
    }

    // Create geometries and collision objects
    m_geoms.reserve(m_entities.size());
    m_objects.reserve(m_entities.size());
    m_kinds.reserve(m_entities.size());
    std::vector<coal::CollisionObject*> raw;
    raw.reserve(m_entities.size());
    for (auto e : m_entities) {
        ColliderKind kind = inferKind_(e);
        m_kinds.push_back(kind);
        auto geom = makeGeometryFor_(kind, e);
        m_geoms.push_back(geom);
        auto obj = std::make_unique<coal::CollisionObject>(geom, /*compute_local_aabb*/ true);
        obj->setTransform(makeTfFromEcs(reg, e));
        obj->computeAABB();

        // set stable indices as user data for fast reverse mapping
        const std::size_t idx = m_objects.size(); // index after previous pushes
        geom->setUserData(reinterpret_cast<void*>(static_cast<uintptr_t>(idx + 1))); // +1 to avoid null
        obj->setUserData(reinterpret_cast<void*>(static_cast<uintptr_t>(idx + 1)));

        raw.push_back(obj.get());
        m_objects.push_back(std::move(obj));
    }

    // Register in broadphase and setup
    m_broadphase->registerObjects(raw);
    m_broadphase->setup();
}

void CollisionCoal::applyTransforms() {
    if (!m_sys) return;
    // Lazily build scene on first use or after clear
    if (!m_broadphase || m_objects.empty()) {
        rebuild();
        if (!m_broadphase) return;
    }
    std::vector<coal::CollisionObject*> updated; updated.reserve(m_objects.size());
    for (std::size_t i = 0; i < m_objects.size(); ++i) {
        entt::entity e = m_entities[i];
        auto* obj = m_objects[i].get();
        // Do NOT apply object transform to halfspaces; their plane is encoded directly in the geometry.
        if (m_kinds[i] == ColliderKind::Halfspace) {
            coal::Transform3s X; X.setIdentity(); obj->setTransform(X);
        } else {
            obj->setTransform(makeTfFromEcs(m_sys->ecs(), e));
        }
        obj->computeAABB();
        updated.push_back(obj);
    }
    m_broadphase->update(updated);
}

std::vector<Contact> CollisionCoal::detectAll() const {
    ContactMap mapCurr;
    if (!m_sys) return {};
    // Ensure the scene exists on first use
    if (!m_broadphase || m_objects.empty()) {
        auto* self = const_cast<CollisionCoal*>(this);
        self->rebuild();
    }
    if (!m_broadphase) return {};

    // First collect potentially colliding pairs from broadphase.
    // We then run narrowphase per pair to get a pair-specific CollisionResult and expand to contact patches.
    coal::CollisionCallBackCollect collect_cb(/*max_size*/ m_objects.size() * m_objects.size());
    collect_cb.init();
    m_broadphase->collide(&collect_cb);

    const auto& reg = m_sys->ecs();
    const std::string frictionCombine = m_sys->config().friction_combine;
    const auto& pairs = collect_cb.getCollisionPairs();
    // Reserve buckets heuristically
    mapCurr.reserve(pairs.size());

    // Prepare reusable request/result for collide and patch computation
    coal::CollisionRequest creq;
    creq.num_max_contacts = 1024;      // per pair cap
    creq.security_margin = -1e-6;      // include touching
    creq.enable_contact = true;
    coal::CollisionResult cres;

    // Contact patch request/result
    coal::ContactPatchRequest patch_req(/*max_num_patch*/ 4);
    coal::ContactPatchResult patch_res(patch_req);

    for (const auto& pr : pairs) {
        auto* o1 = pr.first;
        auto* o2 = pr.second;
        if (!o1 || !o2) continue;

        // Run narrowphase for this pair
        cres.clear();
        coal::collide(o1, o2, creq, cres);
        if (cres.numContacts() == 0) continue;

        // Identify indices/entities back from user data
        auto idx1p = reinterpret_cast<std::uintptr_t>(o1->getUserData());
        auto idx2p = reinterpret_cast<std::uintptr_t>(o2->getUserData());
        if (!(idx1p > 0 && (idx1p - 1) < m_entities.size() && idx2p > 0 && (idx2p - 1) < m_entities.size())) continue;
        const std::size_t iA = idx1p - 1, iB = idx2p - 1;
        entt::entity ea = m_entities[iA];
        entt::entity eb = m_entities[iB];

        // Compute contact patches for this pair based on the pair CollisionResult,
        // then append either expanded patch vertices or raw contacts via a single helper.
        patch_res.clear();
        coal::computeContactPatch(o1, o2, cres, patch_req, patch_res);
        appendContactsFromPair(reg, ea, eb, cres, patch_res, mapCurr, frictionCombine);
    }
    // Try to match current contacts to previous generation per pair (for warmstarting)
    const real_t maxMatchDist = m_sys ? m_sys->config().collision_match_max_dist : (real_t)0.02;
    std::size_t totalContacts = 0;
    std::size_t totalMatched = 0;
    for (auto& kv : mapCurr) {
        auto itPrev = m_prevContactMap.find(kv.first);
        if (itPrev != m_prevContactMap.end()) {
            matchContactsForPair(itPrev->second, kv.second, maxMatchDist);
        } else {
            for (auto& c : kv.second) c.prev_global_out_index = -1;
        }
        // Accumulate simple stats
        totalContacts += kv.second.size();
        for (const auto& c : kv.second) if (c.prev_global_out_index >= 0) ++totalMatched;
    }

    // Flatten map to vector (stable-ish order by map iteration) and assign global indices
    std::vector<Contact> out; out.reserve(1024);
    for (auto& kv : mapCurr) {
        auto& list = kv.second;
        for (std::size_t i = 0; i < list.size(); ++i) {
            list[i].global_out_index = (int)out.size();
            out.emplace_back(list[i]);
        }
    }
    // Optional simple diagnostics
    if (m_sys && m_sys->config().debug_rb) {
        std::printf("[Collision] matched %zu/%zu contacts (threshold = %.4g)\n",
                    totalMatched, totalContacts, (double)maxMatchDist);
    }

    // Update previous contact map for potential warmstarting (next iteration)
    m_prevContactMap = std::move(mapCurr);
    return out;
}

void CollisionCoal::clear() {
    if (m_broadphase) m_broadphase->clear();
    m_objects.clear();
    m_geoms.clear();
    m_entities.clear();
    m_index_from_entity.clear();
}

} // namespace cardillo::collision
