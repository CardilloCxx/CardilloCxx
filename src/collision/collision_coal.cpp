// COAL-based collision manager implementation (persistent, broadphase-backed)
#include "collision_coal.hpp"
#include "../config/config.hpp"
#include "../misc/timings/TimingManager.hpp"
#include "../physics/world.hpp"
#include "../rigid_body/transformations.hpp"
#include "types.hpp"

#include <coal/broadphase/broadphase.h>
#include <coal/broadphase/default_broadphase_callbacks.h>
#include <coal/collision.h>
#include <coal/collision_data.h>
#include <coal/contact_patch.h>
#include <coal/math/transform.h>
#include <coal/shape/geometric_shapes.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include "tangent_frame.hpp"

namespace cardillo::collision {

CollisionCoal::CollisionCoal(cardillo::World& world, cardillo::misc::TimingManager* timings, cardillo::config::Config& cfg)
    : m_world(&world), m_timings(timings), m_cfg(cfg), m_contactTracker(cfg, timings) {
    ensureBroadphaseFromConfig_();
}
CollisionCoal::~CollisionCoal() = default;

namespace {
inline coal::Quatf toCoalQuat(const Quaternion4r& q) {
    return coal::Quatf(q.w(), q.x(), q.y(), q.z());
}
inline coal::Vec3s toCoalVec3(const Vector3r& v) {
    return coal::Vec3s(v.x(), v.y(), v.z());
}
inline coal::Transform3s makeTfFromEcs(const entt::registry& reg, entt::entity e) {
    coal::Transform3s X;
    X.setIdentity();
    Vector3r x = Vector3r::Zero();
    Quaternion4r q = Quaternion4r::Identity();
    if (reg.any_of<cardillo::C_Position3>(e)) x = reg.get<cardillo::C_Position3>(e).value;
    if (reg.any_of<cardillo::C_Orientation>(e)) q = reg.get<cardillo::C_Orientation>(e).value;
    if (reg.any_of<cardillo::C_RB_Cube>(e)) {
        const auto& cb = reg.get<cardillo::C_RB_Cube>(e);
        // Apply local cube orientation on top of body orientation
        q = q * cb.q;
        x += q * cb.center;  // center expressed in cube-local frame
    } else if (reg.any_of<cardillo::C_Cube>(e)) {
        // If only visual cube exists, still honor its center/q
        const auto& cb = reg.get<cardillo::C_Cube>(e);
        q = q * cb.q;
        x += q * cb.center;
    }
    q.normalize();
    X.setTranslation(toCoalVec3(x));
    X.setQuatRotation(toCoalQuat(q));
    return X;
}

// Deduplicate contacts within a pair by proximity.
// Keeps the deepest-penetrating contacts and drops others closer than minDist to any kept one.
static inline void dedupeContactsForPair(ContactList& list, real_t minDist) {
    if (list.size() <= 1 || minDist <= (real_t)0) return;
    const real_t min2 = minDist * minDist;
    // Sort indices by penetration descending (keep deeper contacts first)
    std::vector<std::size_t> idx(list.size());
    for (std::size_t i = 0; i < list.size(); ++i) idx[i] = i;
    std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) { return list[a].penetration > list[b].penetration; });
    std::vector<Contact> kept;
    kept.reserve(list.size());
    for (std::size_t k = 0; k < idx.size(); ++k) {
        const Contact& c = list[idx[k]];
        bool tooClose = false;
        for (const auto& sel : kept) {
            const Vector3r d = c.point - sel.point;
            if (d.squaredNorm() <= min2) {
                tooClose = true;
                break;
            }
        }
        if (!tooClose) kept.emplace_back(c);
    }
    list.swap(kept);
}

// Helper to populate a Contact from world data and entities
static inline real_t combineFrictionMu(const entt::registry& reg, entt::entity a, entt::entity b, const std::string& method) {
    auto getMu = [&](entt::entity e) -> real_t {
        if (reg.any_of<cardillo::C_Friction>(e)) {
            return std::max<real_t>((real_t)0, reg.get<cardillo::C_Friction>(e).mu);
        }
        return (real_t)0;
    };
    const real_t muA = getMu(a);
    const real_t muB = getMu(b);
    if (method == "arith") return (muA + muB) * (real_t)0.5;
    if (method == "geom") return (real_t)std::sqrt(std::max<real_t>((real_t)0, muA * muB));
    return std::min(muA, muB);  // default: min
}

// Helper to populate a Contact from world data and entities. Takes the bodies' RigidState and
// combined friction coefficient as inputs rather than recomputing them -- both depend only on the
// (ea, eb) pair, not on the specific contact point, so callers compute them once per pair and
// reuse across every point in a multi-point manifold (see appendContactsFromPair below).
inline Contact makeContact(entt::entity ea, entt::entity eb, const cardillo::RigidBody::RigidState& stateA, const cardillo::RigidBody::RigidState& stateB, real_t friction_mu, Vector3r p1W,
                            Vector3r p2W, Vector3r nN, real_t depth) {
    Contact c{};
    c.a = ea;
    c.b = eb;

    // Same class of bug as the normal below, hitting the nearest/patch points instead: COAL can
    // return a non-finite point for one (or both) sides of an exactly-touching/coincident pair.
    // Sanitize before anything derives from them (c.point, and from there pointA_body/pointB_body
    // -- which feed the angular row of buildContactRowByDof() for any 6-DOF body, so a NaN point
    // silently poisons that body's Jacobian even when the normal itself is fine).
    if (!p1W.allFinite()) p1W = p2W.allFinite() ? p2W : stateA.position;
    if (!p2W.allFinite()) p2W = p1W.allFinite() ? p1W : stateB.position;

    // Normalize normal and construct stable tangent frame. `nlen > 0` alone is not a sufficient
    // guard: for two primitives that are exactly touching or exactly coincident (e.g. a chain of
    // abutting spheres spaced at exactly the sum of their radii -- see the hangbridge scene's
    // closely-packed rope nodes), COAL's own narrowphase can return an already-NaN normal, and
    // `NaN > 0` is false, so that case fell through this check untouched and propagated NaN into
    // the whole contact (and from there into the solver's rhs/Delassus blocks). Explicitly check
    // finiteness and fall back, in order, to the point-separation direction, then the direction
    // between body centers, then an arbitrary axis -- never leave a non-finite normal in `c`.
    real_t nlen = nN.norm();
    if (nlen > (real_t)0 && std::isfinite(nlen)) {
        nN /= nlen;
    } else {
        Vector3r fallback = p2W - p1W;
        real_t flen = fallback.norm();
        if (!(flen > (real_t)0) || !std::isfinite(flen)) {
            fallback = stateB.position - stateA.position;
            flen = fallback.norm();
        }
        nN = (flen > (real_t)0 && std::isfinite(flen)) ? (fallback / flen) : Vector3r(0, 0, 1);
    }

    c.normal = nN;  // convention: from A -> B
    cardillo::collision::tangentFrameFromNormal(nN, c.tangent1, c.tangent2);
    c.point = (p1W + p2W) * (real_t)0.5;
    c.penetration = std::max<real_t>(0.0, depth);
    const cardillo::RigidBody::RigidState inertial = cardillo::RigidBody::RigidState::inertial();
    c.pointA_body = cardillo::transform::point(c.point, inertial, stateA);
    c.normalA_body = cardillo::transform::direction(c.normal, inertial, stateA);
    c.pointB_body = cardillo::transform::point(c.point, inertial, stateB);
    c.normalB_body = cardillo::transform::direction(c.normal, inertial, stateB);
    // Tangents in body frames
    c.tangent1A_body = cardillo::transform::direction(c.tangent1, inertial, stateA);
    c.tangent2A_body = cardillo::transform::direction(c.tangent2, inertial, stateA);
    c.tangent1B_body = cardillo::transform::direction(c.tangent1, inertial, stateB);
    c.tangent2B_body = cardillo::transform::direction(c.tangent2, inertial, stateB);
    c.friction_mu = friction_mu;
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

// Append contact-manifold geometry for a pair into `manifolds` -- mirrors appendContactsFromPair's
// two branches (patch expansion vs. raw fallback) below so this stream is always at least as
// complete as the flattened Contact list: whenever patch computation yields nothing for a pair
// (e.g. a too-degenerate vertex-vertex contact), appendContactsFromPair falls back to the raw
// narrowphase contact(s), and this must fall back the same way or manifolds would silently miss
// pairs that the flattened list still shows. Unlike appendContactsFromPair, points stay grouped
// per patch (a manifold can have several points) instead of being flattened into independent
// per-vertex Contacts. Reuses makeContact() (below) for its point/normal/tangent-frame
// sanitization (NaN/degenerate guards) instead of re-deriving that logic here.
inline void appendManifoldsForPair(entt::registry& reg, entt::entity ea, entt::entity eb, const coal::CollisionResult& cres, const coal::ContactPatchResult& patch_res,
                                    std::vector<ContactManifold>& manifolds, const std::string& frictionCombine, bool usePatchVertices);

// Append contacts for a pair into a map: prefer patch expansion; fall back to raw contacts
inline std::size_t appendContactsFromPair(entt::registry& reg, entt::entity ea, entt::entity eb, const coal::CollisionResult& cres, const coal::ContactPatchResult& patch_res, ContactMap& outMap,
                                          const std::string& frictionCombine, bool usePatchVertices) {
    // Depend only on the (ea, eb) pair, not on the individual contact point -- compute once per
    // pair instead of once per point (a manifold can have several points per pair).
    const cardillo::RigidBody::RigidState stateA = cardillo::RigidBody::getState(reg, ea);
    const cardillo::RigidBody::RigidState stateB = cardillo::RigidBody::getState(reg, eb);
    const real_t friction_mu = combineFrictionMu(reg, ea, eb, frictionCombine);

    std::size_t appended = 0;
    if (usePatchVertices && patch_res.numContactPatches() > 0) {
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
                addContactToMap(outMap, makeContact(ea, eb, stateA, stateB, friction_mu, p1W, p2W, nW, depth));
                ++appended;
            }
        }
        return appended;
    }

    // Fallback: emit raw contacts from narrowphase result
    for (int k = 0; k < cres.numContacts(); ++k) {
        const coal::Contact c0 = cres.getContact(k);
        const Vector3r nW(c0.normal.x(), c0.normal.y(), c0.normal.z());
        const Vector3r p1W(c0.nearest_points[0].x(), c0.nearest_points[0].y(), c0.nearest_points[0].z());
        const Vector3r p2W(c0.nearest_points[1].x(), c0.nearest_points[1].y(), c0.nearest_points[1].z());
        addContactToMap(outMap, makeContact(ea, eb, stateA, stateB, friction_mu, p1W, p2W, nW, (real_t)c0.penetration_depth));
        ++appended;
    }
    return appended;
}

inline void appendManifoldsForPair(entt::registry& reg, entt::entity ea, entt::entity eb, const coal::CollisionResult& cres, const coal::ContactPatchResult& patch_res,
                                    std::vector<ContactManifold>& manifolds, const std::string& frictionCombine, bool usePatchVertices) {
    const cardillo::RigidBody::RigidState stateA = cardillo::RigidBody::getState(reg, ea);
    const cardillo::RigidBody::RigidState stateB = cardillo::RigidBody::getState(reg, eb);
    const real_t friction_mu = combineFrictionMu(reg, ea, eb, frictionCombine);

    if (usePatchVertices && patch_res.numContactPatches() > 0) {
        for (std::size_t ip = 0; ip < patch_res.numContactPatches(); ++ip) {
            const coal::ContactPatch& patch = patch_res.getContactPatch(ip);
            const std::size_t m = patch.size();
            if (m == 0) continue;
            const Vector3r nW(patch.getNormal().x(), patch.getNormal().y(), patch.getNormal().z());
            const real_t depth = (real_t)patch.penetration_depth;

            ContactManifold cm;
            cm.a = ea;
            cm.b = eb;
            cm.friction_mu = friction_mu;
            cm.points.reserve(m);
            for (std::size_t iv = 0; iv < m; ++iv) {
                const auto p1Wc = patch.getPointShape1(iv);
                const auto p2Wc = patch.getPointShape2(iv);
                const Vector3r p1W(p1Wc.x(), p1Wc.y(), p1Wc.z());
                const Vector3r p2W(p2Wc.x(), p2Wc.y(), p2Wc.z());
                const Contact c = makeContact(ea, eb, stateA, stateB, friction_mu, p1W, p2W, nW, depth);
                if (iv == 0) {
                    cm.normal = c.normal;
                    cm.tangent1 = c.tangent1;
                    cm.tangent2 = c.tangent2;
                    cm.penetration = c.penetration;
                }
                ContactManifold::Point pt;
                pt.position = c.point;
                cm.points.push_back(pt);
            }
            manifolds.push_back(std::move(cm));
        }
        return;
    }

    // Fallback: same raw narrowphase contacts appendContactsFromPair uses when patch computation
    // yielded nothing for this pair -- one single-point manifold per raw contact.
    for (int k = 0; k < cres.numContacts(); ++k) {
        const coal::Contact c0 = cres.getContact(k);
        const Vector3r nW(c0.normal.x(), c0.normal.y(), c0.normal.z());
        const Vector3r p1W(c0.nearest_points[0].x(), c0.nearest_points[0].y(), c0.nearest_points[0].z());
        const Vector3r p2W(c0.nearest_points[1].x(), c0.nearest_points[1].y(), c0.nearest_points[1].z());
        const Contact c = makeContact(ea, eb, stateA, stateB, friction_mu, p1W, p2W, nW, (real_t)c0.penetration_depth);

        ContactManifold cm;
        cm.a = ea;
        cm.b = eb;
        cm.normal = c.normal;
        cm.tangent1 = c.tangent1;
        cm.tangent2 = c.tangent2;
        cm.penetration = c.penetration;
        cm.friction_mu = friction_mu;
        ContactManifold::Point pt;
        pt.position = c.point;
        cm.points.push_back(pt);
        manifolds.push_back(std::move(cm));
    }
}
}  // namespace

void CollisionCoal::ensureBroadphaseFromConfig_() {
    if (m_broadphase) return;
    // Select based on config
    std::string bp = m_world ? m_cfg.collision_broadphase : std::string("dynamic_aabb");
    auto toLower = [](std::string s) {
        for (auto& c : s) c = (char)std::tolower((unsigned char)c);
        return s;
    };
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
    const auto& reg = m_world->ecs();
    if (reg.any_of<cardillo::C_RB_Cube>(e)) return ColliderKind::Box;
    if (reg.any_of<cardillo::C_RB_Capsule>(e)) return ColliderKind::Capsule;
    if (reg.any_of<cardillo::C_RB_Cylinder>(e)) return ColliderKind::Cylinder;
    if (reg.any_of<cardillo::C_RB_Cone>(e)) return ColliderKind::Cone;
    if (reg.any_of<cardillo::C_RB_Plane>(e)) return ColliderKind::Halfspace;
    if ((reg.any_of<cardillo::C_PointMassTag>(e) || reg.any_of<cardillo::C_RB_Sphere>(e)) && reg.any_of<cardillo::C_Radius>(e)) return ColliderKind::Sphere;
    if (reg.any_of<cardillo::C_RB_Mesh>(e) && reg.any_of<cardillo::C_Mesh>(e)) return ColliderKind::Mesh;

    throw std::runtime_error("CollisionCoal: unsupported collider entity; add appropriate tag/geometry.");
}

std::shared_ptr<coal::CollisionGeometry> CollisionCoal::makeGeometryFor_(ColliderKind kind, entt::entity e) const {
    const auto& reg = m_world->ecs();
    switch (kind) {
        case ColliderKind::Box: {
            const auto& he = reg.get<cardillo::C_RB_Cube>(e).halfExtents;
            return std::make_shared<coal::Box>(he.x() * 2.0, he.y() * 2.0, he.z() * 2.0);
        }
        case ColliderKind::Halfspace: {
            const auto& plane = reg.get<cardillo::C_RB_Plane>(e);
            Vector3r n = plane.normal.normalized();
            const auto& x = reg.get<cardillo::C_Position3>(e).value;
            real_t d = n.dot(x);
            return std::make_shared<coal::Halfspace>(toCoalVec3(n), (coal::CoalScalar)d);
        }
        case ColliderKind::Sphere: {
            const auto r = reg.get<cardillo::C_Radius>(e).r;
            return std::make_shared<coal::Sphere>((coal::CoalScalar)r);
        }
        case ColliderKind::Capsule: {
            const auto& cap = reg.get<cardillo::C_RB_Capsule>(e);
            return std::make_shared<coal::Capsule>((coal::CoalScalar)cap.radius, (coal::CoalScalar)(cap.halfLength * 2));
        }
        case ColliderKind::Cylinder: {
            const auto& cyl = reg.get<cardillo::C_RB_Cylinder>(e);
            return std::make_shared<coal::Cylinder>((coal::CoalScalar)cyl.radius, (coal::CoalScalar)(cyl.halfLength * 2));
        }
        case ColliderKind::Cone: {
            const auto& cone = reg.get<cardillo::C_RB_Cone>(e);
            return std::make_shared<coal::Cone>((coal::CoalScalar)cone.radius, (coal::CoalScalar)cone.height);
        }
        case ColliderKind::Mesh: {
            const auto& asset = m_world->getMeshAsset(e);
            if (!asset.bvh) throw std::runtime_error("COAL MeshLoader failed to load BVH for mesh entity");
            return asset.bvh;
        }
    }
    throw std::runtime_error("CollisionCoal: unknown collider kind");
}

void CollisionCoal::rebuild() {
    if (!m_world) return;
    auto sc = m_timings->scope(cardillo::misc::TimingManager::TimerId::CollisionBroadphase);
    clear();
    ensureBroadphaseFromConfig_();
    auto& reg = m_world->ecs();

    // Gather all collidable entities in a single pass
    m_entities.clear();
    m_index_from_entity.clear();
    {
        auto view = reg.view<cardillo::C_Collidable, cardillo::C_Position3>();
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
        const std::size_t idx = m_objects.size();                                     // index after previous pushes
        geom->setUserData(reinterpret_cast<void*>(static_cast<uintptr_t>(idx + 1)));  // +1 to avoid null
        obj->setUserData(reinterpret_cast<void*>(static_cast<uintptr_t>(idx + 1)));

        raw.push_back(obj.get());
        m_objects.push_back(std::move(obj));
    }

    // Register in broadphase and setup
    m_broadphase->registerObjects(raw);
    m_broadphase->setup();
}

void CollisionCoal::applyTransforms() {
    if (!m_world) return;
    auto sc = m_timings->scope(cardillo::misc::TimingManager::TimerId::CollisionBroadphase);
    // Lazily build scene on first use or after clear
    if (!m_broadphase || m_objects.empty()) {
        rebuild();
        if (!m_broadphase) return;
    }
    std::vector<coal::CollisionObject*> updated;
    updated.reserve(m_objects.size());
    for (std::size_t i = 0; i < m_objects.size(); ++i) {
        entt::entity e = m_entities[i];
        auto* obj = m_objects[i].get();
        const bool isDynamic = m_world->ecs().any_of<cardillo::C_PhysicsObject, cardillo::C_StaticTrajectory>(e);
        if (!isDynamic) {
            // Static objects: transform was set at rebuild; no need to recompute each step
            continue;
        }
        // Do NOT apply object transform to halfspaces; their plane is encoded directly in the
        // geometry.
        if (m_kinds[i] == ColliderKind::Halfspace) {
            coal::Transform3s X;
            X.setIdentity();
            obj->setTransform(X);
        } else {
            obj->setTransform(makeTfFromEcs(m_world->ecs(), e));
        }
        obj->computeAABB();
        updated.push_back(obj);
    }
    m_broadphase->update(updated);
}

std::vector<Contact>& CollisionCoal::detectAll() {
    ContactMap& mapCurr = m_mapCurr;
    mapCurr.clear();
    m_contactManifolds.clear();
    if (!m_world) {
        m_prev_flattened.clear();
        m_flattened.clear();
        return m_flattened;
    }
    // Ensure the scene exists on first use
    if (!m_broadphase || m_objects.empty()) {
        rebuild();
    }

    auto& reg = m_world->ecs();

    std::vector<coal::CollisionCallBackCollect::CollisionPair>& pairs = m_pairs;
    pairs.clear();
    {
        auto sc2 = m_timings->scope(cardillo::misc::TimingManager::TimerId::CollisionBroadphase);
        const std::size_t nObj = m_objects.size();
        const std::size_t cbReserve = std::max<std::size_t>(nObj * 16, 256);
        coal::CollisionCallBackCollect collect_cb(/*max_size*/ cbReserve);
        collect_cb.init();
        m_broadphase->collide(&collect_cb);
        pairs = collect_cb.getCollisionPairs();
    }

    // Reserve buckets heuristically
    mapCurr.reserve(pairs.size());

    if (pairs.empty()) {
        if (m_world && m_cfg.debug_rb) {
            std::printf("[Collision] no potential collision pairs found in broadphase.\n");
        }
        // clear previous/current buffers
        m_prev_flattened.clear();
        m_flattened.clear();
        return m_flattened;
    }

    // Prepare reusable request/result for collide and patch computation
    coal::CollisionRequest creq;
    creq.num_max_contacts = m_world ? m_cfg.collision_max_raw_contacts : 1024;  // per pair cap
    // Use configurable security margin (>=0). A small positive value improves robustness for thin
    // triangles
    creq.security_margin = m_world ? m_cfg.collision_security_margin : (real_t)1e-5;
    creq.enable_contact = true;
    coal::CollisionResult cres;

    // Contact patch request/result
    const std::size_t max_patch_req = m_world ? m_cfg.collision_max_patches : (std::size_t)4;
    coal::ContactPatchRequest patch_req(/*max_num_patch*/ max_patch_req);
    coal::ContactPatchResult patch_res(patch_req);
    const bool usePatchVertices = m_world ? m_cfg.collision_use_patch_vertices : true;

    const std::string frictionCombine = m_cfg.friction_combine;

    for (const auto& pr : pairs) {
        auto sc_n = m_timings->scope(cardillo::misc::TimingManager::TimerId::CollisionNarrowphase);
        auto* o1 = pr.first;
        auto* o2 = pr.second;
        if (!o1 || !o2) continue;

        // Run narrowphase for this pair
        entt::entity ea, eb;
        {
            auto sc_c = m_timings->scope(cardillo::misc::TimingManager::TimerId::CollisionNarrowphaseCollide);

            // Identify indices/entities back from user data
            auto idx1p = reinterpret_cast<std::uintptr_t>(o1->getUserData());
            auto idx2p = reinterpret_cast<std::uintptr_t>(o2->getUserData());
            if (!(idx1p > 0 && (idx1p - 1) < m_entities.size() && idx2p > 0 && (idx2p - 1) < m_entities.size())) continue;
            std::size_t iA = idx1p - 1, iB = idx2p - 1;

            if (iB < iA) {
                std::swap(o1, o2);
                std::swap(iA, iB);
            }

            ea = m_entities[iA];
            eb = m_entities[iB];

            if (RigidBody::isStatic(reg, ea) && RigidBody::isStatic(reg, eb)) continue;
            if (isPairDisabled(ea, eb)) continue;

            cres.clear();
            coal::collide(o1, o2, creq, cres);
            if (cres.numContacts() == 0) continue;
        }

        // Compute contact patches for this pair
        {
            auto sc_prep = m_timings->scope(cardillo::misc::TimingManager::TimerId::CollisionMakeContactPatch);
            patch_res.clear();
            coal::computeContactPatch(o1, o2, cres, patch_req, patch_res);
        }

        // Append contacts (flattened, for the solver) and manifolds (grouped per patch, for
        // export/visualization) from this pair. appendManifoldsForPair mirrors
        // appendContactsFromPair's patch/fallback branches so the manifold stream is always at
        // least as complete as the flattened contact list -- see its doc comment above.
        {
            auto sc_p = m_timings->scope(cardillo::misc::TimingManager::TimerId::CollisionMakeContact);
            appendContactsFromPair(reg, ea, eb, cres, patch_res, mapCurr, frictionCombine, usePatchVertices);
            appendManifoldsForPair(reg, ea, eb, cres, patch_res, m_contactManifolds, frictionCombine, usePatchVertices);
        }
    }

    // Dedupe contacts within each pair by proximity
    const real_t minPairDist = m_world ? m_cfg.collision_min_pair_contact_distance : (real_t)0.0;
    if (minPairDist > (real_t)0) {
        for (auto& kv : mapCurr) {
            dedupeContactsForPair(kv.second, minPairDist);
        }
    }

    // Move previous authoritative flattened buffer aside so tracker can use it
    m_prev_flattened = std::move(m_flattened);
    m_flattened.clear();

    // Copy last_impulse values from authoritative previous flattened buffer into current ContactMap.
    // Reserve based on the previous step's actual contact count rather than a fixed guess, since
    // contact counts are typically stable step-to-step outside of impact events.
    m_flattened.reserve(std::max<std::size_t>(m_prev_flattened.size(), 64));
    int globalIndex = 0;
    for (auto& [key, list] : mapCurr) {
        for (int i = 0; i < list.size(); ++i) {
            Contact& c = list[i];
            m_flattened.emplace_back(std::move(c));
            m_flattened.back().global_out_index = globalIndex++;
            mapCurr[key][i].global_out_index = m_flattened.back().global_out_index;
        }
    }

    if (m_cfg.pj_warmstart) {
        m_contactTracker.registerNextContacts(mapCurr);
        m_contactTracker.applyPrevImpulses(m_flattened, m_prev_flattened, mapCurr);
    }

    return m_flattened;
}

void CollisionCoal::clear() {
    if (m_broadphase) m_broadphase->clear();
    m_objects.clear();
    m_geoms.clear();
    m_entities.clear();
    m_index_from_entity.clear();
}

void CollisionCoal::disablePair(entt::entity a, entt::entity b) {
    m_disabledPairs.insert(ContactPairKey::make(a, b));
}

void CollisionCoal::enablePair(entt::entity a, entt::entity b) {
    m_disabledPairs.erase(ContactPairKey::make(a, b));
}

bool CollisionCoal::isPairDisabled(entt::entity a, entt::entity b) const {
    auto sc3 = m_timings->scope(cardillo::misc::TimingManager::TimerId::DisableCollisionPairs);
    return m_disabledPairs.find(ContactPairKey::make(a, b)) != m_disabledPairs.end();
}

}  // namespace cardillo::collision
