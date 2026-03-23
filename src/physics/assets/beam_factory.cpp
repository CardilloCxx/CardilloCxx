#include "beam_factory.hpp"

#include "../constraints/constraints.hpp"

#include <cmath>

namespace cardillo::physics {
namespace {

struct BeamSample {
    Vector3r position;
    Vector3r tangent;
    Vector3r normal;
    Vector3r binormal;
    real_t segLen;
};

Matrix33r makeFrameFromTangentLocal(const Vector3r& tangent) {
    Vector3r T = tangent.normalized();
    Vector3r up = std::abs(T.z()) < (real_t)0.9 ? Vector3r(0, 0, 1) : Vector3r(0, 1, 0);
    Vector3r N = up.cross(T).normalized();
    if (N.squaredNorm() < (real_t)1e-12) {
        up = Vector3r(1, 0, 0);
        N = up.cross(T).normalized();
    }
    Vector3r B = T.cross(N).normalized();
    Matrix33r M;
    M.col(0) = T;
    M.col(1) = N;
    M.col(2) = B;
    return M;
}

std::pair<entt::entity, entt::entity> buildBeamFromSamples(
    World& sys,
    const std::vector<BeamSample>& samples,
    bool loop,
    const BeamCrossSection& section,
    const BeamSpringParams& springs,
    const RigidState& stateDefaults,
    const RigidProps& propsDefaults,
    const Vector3r& splineCOMWorld) {
    if (samples.empty()) return {entt::null, entt::null};

    real_t totalLen = (real_t)0;
    for (const auto& s : samples) totalLen += s.segLen;
    if (totalLen <= (real_t)0) totalLen = (real_t)1;

    Matrix33r Rshape = Matrix33r::Identity();
    if (section.type == BeamBodyType::Capsule ||
        section.type == BeamBodyType::Cylinder) {
        Rshape = Quaternion4r::FromTwoVectors(Vector3r::UnitZ(), Vector3r::UnitX()).toRotationMatrix();
    }

    entt::entity root = entt::null;
    entt::entity prev = entt::null;
    entt::entity end = entt::null;

    const Matrix33r Rbody = stateDefaults.orientation.toRotationMatrix();
    const Vector3r worldCOM = splineCOMWorld + stateDefaults.position;
    const Vector3r v_body_world = Rbody * stateDefaults.linearVelocity;
    const Vector3r w_body_world = Rbody * stateDefaults.angularVelocity;

    Quaternion4r q_prev = Quaternion4r::Identity();

    for (const auto& s : samples) {
        const real_t segLen = s.segLen;

        RigidShape shape;
        if (section.type == BeamBodyType::Cube) {
            shape = CubeShape(
                Vector3r(segLen * (real_t)0.5, section.width * (real_t)0.5, section.height * (real_t)0.5));
        } else if (section.type == BeamBodyType::Cylinder) {
            const real_t r = std::min(section.width, section.height) * (real_t)0.5;
            shape = CylinderShape(r, segLen * (real_t)0.5);
        } else {
            const real_t r = std::min(section.width, section.height) * (real_t)0.5;
            shape = CapsuleShape(r, segLen * (real_t)0.5);
        }

        RigidProps segProps = propsDefaults;
        real_t massPerSegment = (real_t)0;
        if (propsDefaults.mass.has_value()) {
            massPerSegment = *propsDefaults.mass * (segLen / totalLen);
        } else if (propsDefaults.density.has_value()) {
            massPerSegment = *propsDefaults.density * (section.area() * segLen);
        }
        segProps.mass = (massPerSegment > (real_t)0) ? std::optional<real_t>(massPerSegment) : std::nullopt;

        Matrix33r Rlocal;
        if (s.normal.squaredNorm() > (real_t)0 && s.binormal.squaredNorm() > (real_t)0) {
            Rlocal.col(0) = s.tangent.normalized();
            Rlocal.col(1) = s.normal.normalized();
            Rlocal.col(2) = s.binormal.normalized();
        } else {
            Rlocal = makeFrameFromTangentLocal(s.tangent);
        }

        const Matrix33r Rworld = Rbody * Rlocal * Rshape;
        Quaternion4r qworld(Rworld);
        qworld.normalize();

        qworld = World::alignQuaternionTo(qworld, q_prev);
        q_prev = qworld;

        const Vector3r worldPos = splineCOMWorld + stateDefaults.position + Rbody * (s.position - splineCOMWorld);
        const Vector3r v_world = v_body_world + w_body_world.cross(worldPos - worldCOM);

        RigidState segState;
        segState.position = worldPos;
        segState.orientation = qworld;
        segState.linearVelocity = v_world;
        segState.angularVelocity = Rlocal.transpose() * stateDefaults.angularVelocity;

        const entt::entity cur = sys.addRigidBody(shape, segState, segProps);

        TopologyComponents::BeamElement be_cur;
        be_cur.l0 = segLen;
        be_cur.l = segLen;
        if (prev != entt::null) {
            be_cur.prev = prev;
            if (!sys.ecs().any_of<TopologyComponents::BeamElement>(prev)) {
                TopologyComponents::BeamElement be_prev;
                be_prev.l0 = segLen;
                be_prev.l = segLen;
                be_prev.next = cur;
                sys.ecs().emplace<TopologyComponents::BeamElement>(prev, be_prev);
            } else {
                sys.ecs().get<TopologyComponents::BeamElement>(prev).next = cur;
            }
        }
        sys.ecs().emplace<TopologyComponents::BeamElement>(cur, be_cur);

        if (prev != entt::null) {
            sys.addBeamConstraint(prev, cur, springs, section);
            sys.disableCollisionBetween(prev, cur);
        }

        if (root == entt::null) root = cur;
        prev = cur;
        end = cur;
    }

    if (loop && root != entt::null && end != entt::null && end != root) {
        sys.addBeamConstraint(end, root, springs, section);
        sys.disableCollisionBetween(end, root);
        if (sys.ecs().any_of<TopologyComponents::BeamElement>(end)) {
            sys.ecs().get<TopologyComponents::BeamElement>(end).next = root;
        }
        if (sys.ecs().any_of<TopologyComponents::BeamElement>(root)) {
            sys.ecs().get<TopologyComponents::BeamElement>(root).prev = end;
        }
    }

    return {root, end};
}

} // namespace

std::pair<entt::entity, entt::entity> BeamFactory::createBeam(
    World& system,
    const misc::SplinePattern& spline,
    const BeamCrossSection& section,
    const BeamSpringParams& springs,
    const RigidState& stateDefaults,
    const RigidProps& propsDefaults,
    size_t segments) {
    const real_t totalLen = spline.totalLength();
    const real_t segLen = totalLen / (real_t)segments;
    const bool endsOnSpline = true;

    std::vector<BeamSample> samples;
    samples.reserve(segments);

    if (endsOnSpline) {
        const real_t minSegLen = (real_t)1e-8;
        const size_t segCount = segments;
        for (size_t i = 0; i < segCount; ++i) {
            real_t alpha0 = (real_t)i / (real_t)segments;
            real_t alpha1 = (real_t)(i + 1) / (real_t)segments;
            if (spline.isLoop() && i + 1 == segCount) alpha1 = (real_t)0;

            const misc::SplineSample si0 = spline.sample(alpha0);
            const misc::SplineSample si1 = spline.sample(alpha1);
            const Vector3r midPos = (si0.position + si1.position) * (real_t)0.5;
            const real_t local_segLen = (si1.position - si0.position).norm();
            if (local_segLen <= minSegLen) continue;

            const Vector3r midTangent = (si1.position - si0.position) / local_segLen;
            const Vector3r midNormal = (si0.normal + si1.normal).normalized();
            const Vector3r midBinormal = (si0.binormal + si1.binormal).normalized();
            samples.push_back(BeamSample{midPos, midTangent, midNormal, midBinormal, local_segLen});
        }
    } else {
        for (size_t i = 0; i <= segments; ++i) {
            const real_t alpha = (real_t)i / (real_t)segments;
            const misc::SplineSample si = spline.sample(alpha);
            samples.push_back(BeamSample{si.position, si.tangent, si.normal, si.binormal, segLen});
        }
    }

    const Vector3r splineCOMWorld = spline.centerOfMass();
    return buildBeamFromSamples(
        system,
        samples,
        spline.isLoop(),
        section,
        springs,
        stateDefaults,
        propsDefaults,
        splineCOMWorld);
}

std::pair<entt::entity, entt::entity> BeamFactory::createBeams(
    World& system,
    const std::vector<const misc::SplinePattern*>& splines,
    const BeamCrossSection& section,
    const BeamSpringParams& springs,
    const RigidState& stateDefaults,
    const RigidProps& propsDefaults,
    size_t segments) {
    real_t totalLen = (real_t)0;
    for (const auto* sp : splines) {
        if (sp) totalLen += sp->totalLength();
    }

    entt::entity first = entt::null;
    entt::entity second = entt::null;
    entt::entity prevEnd = entt::null;

    for (size_t i = 0; i < splines.size(); ++i) {
        const auto pair = createBeam(
            system,
            *splines[i],
            section,
            springs,
            stateDefaults,
            propsDefaults,
            (size_t)(segments * (splines[i]->totalLength() / totalLen)));

        if (first == entt::null) first = pair.first;
        if (prevEnd != entt::null && pair.first != entt::null) {
            if (system.ecs().any_of<MotionComponents::Orientation>(prevEnd) &&
                system.ecs().any_of<MotionComponents::Orientation>(pair.first)) {
                auto& qNext = system.ecs().get<MotionComponents::Orientation>(pair.first).value;
                const auto& qPrev = system.ecs().get<MotionComponents::Orientation>(prevEnd).value;
                qNext = World::alignQuaternionTo(qNext, qPrev);
            }
            system.addRigidConstraint(prevEnd, pair.first);
            system.disableCollisionBetween(prevEnd, pair.first);
        }
        prevEnd = pair.second;
    }

    second = prevEnd;
    return {first, second};
}

} // namespace cardillo::physics
