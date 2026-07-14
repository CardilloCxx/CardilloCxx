#include "beam_factory.hpp"

#include "../../rigid_body/rigid_body.hpp"
#include "../../rigid_body/transformations.hpp"
#include "../assets/body_factory.hpp"
#include "../constraints/constraint_factory.hpp"
#include "../constraints/constraints.hpp"

#include "../../collision/collision_coal.hpp"

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

void appendSplineSamples(const misc::SplinePattern& spline, size_t segments, std::vector<BeamSample>& outSamples) {
    if (segments == 0) return;

    const real_t minSegLen = (real_t)1e-8;
    const size_t segCount = segments;

    for (size_t i = 0; i < segCount; ++i) {
        real_t alpha0 = (real_t)i / (real_t)segCount;
        real_t alpha1 = (real_t)(i + 1) / (real_t)segCount;
        if (spline.isLoop() && i + 1 == segCount) alpha1 = (real_t)0;

        const misc::SplineSample si0 = spline.sample(alpha0);
        const misc::SplineSample si1 = spline.sample(alpha1);
        const real_t localSegLen = (si1.position - si0.position).norm();
        if (localSegLen <= minSegLen) continue;

        const Vector3r midPos = (si0.position + si1.position) * (real_t)0.5;
        const Vector3r midTangent = (si1.position - si0.position) / localSegLen;

        Vector3r n0 = si0.normal;
        Vector3r n1 = si1.normal;
        Vector3r b0 = si0.binormal;
        Vector3r b1 = si1.binormal;

        // Detect and correct sudden 180-degree flips in the frame
        if (n0.dot(n1) < (real_t)0.0) {
            n1 = -n1;
            b1 = -b1; // Flip binormal as well to maintain handedness
        }

        // Average the aligned normals
        Vector3r midNormal = n0 + n1;
        
        // Orthonormalize the normal against the new midTangent via Gram-Schmidt
        midNormal -= midTangent * midTangent.dot(midNormal);
        
        if (midNormal.allFinite() && midNormal.squaredNorm() > minSegLen * minSegLen) {
            midNormal.normalize();
        } else {
            // Fallback: If normal calculation completely fails, generate an arbitrary orthogonal vector
            Vector3r up = (std::abs(midTangent.z()) < (real_t)0.999) ? Vector3r::UnitZ() : Vector3r::UnitX();
            midNormal = midTangent.cross(up).normalized();
        }

        // The binormal is strictly the cross product of tangent and normal
        Vector3r midBinormal = midTangent.cross(midNormal).normalized();

        outSamples.push_back(BeamSample{midPos, midTangent, midNormal, midBinormal, localSegLen});
    }
}

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

Matrix33r makeStableFrameFromSample(const BeamSample& s, const Matrix33r* prevFrame) {
    constexpr real_t eps2 = (real_t)1e-16;

    Vector3r T = s.tangent;
    if (!T.allFinite() || T.squaredNorm() <= eps2) {
        if (prevFrame != nullptr) {
            T = prevFrame->col(0);
        } else {
            T = Vector3r::UnitX();
        }
    }
    T.normalize();

    Vector3r N = s.normal;
    if (!N.allFinite() || N.squaredNorm() <= eps2) {
        if (prevFrame != nullptr) {
            N = prevFrame->col(1);
        } else {
            N = Vector3r::Zero();
        }
    }
    N -= T * T.dot(N);

    if (!N.allFinite() || N.squaredNorm() <= eps2) {
        Vector3r Bcand = s.binormal;
        if (!Bcand.allFinite() || Bcand.squaredNorm() <= eps2) {
            if (prevFrame != nullptr) {
                Bcand = prevFrame->col(2);
            } else {
                Bcand = Vector3r::Zero();
            }
        }
        Bcand -= T * T.dot(Bcand);
        if (Bcand.allFinite() && Bcand.squaredNorm() > eps2) {
            Bcand.normalize();
            N = Bcand.cross(T);
        }
    }

    if (!N.allFinite() || N.squaredNorm() <= eps2) {
        return makeFrameFromTangentLocal(T);
    }

    N.normalize();
    Vector3r B = T.cross(N);
    if (!B.allFinite() || B.squaredNorm() <= eps2) {
        return makeFrameFromTangentLocal(T);
    }
    B.normalize();
    N = B.cross(T).normalized();

    // Keep the normal/binormal orientation continuous along the beam.
    if (prevFrame != nullptr && N.dot(prevFrame->col(1)) < (real_t)0) {
        N = -N;
        B = -B;
    }

    Matrix33r M;
    M.col(0) = T;
    M.col(1) = N;
    M.col(2) = B;
    return M;
}

std::pair<entt::entity, entt::entity> buildBeamFromSamples(World& sys, const std::vector<BeamSample>& samples, bool loop, const BeamCrossSection& section, const BeamSpringParams& springs,
                                                           const RigidState& stateDefaults, const RigidProps& propsDefaults, const Vector3r& splineCOMWorld,
                                                           cardillo::collision::CollisionCoal* collision_mgr) {
    if (samples.empty()) return {entt::null, entt::null};

    real_t totalLen = (real_t)0;
    for (const auto& s : samples) totalLen += s.segLen;
    if (totalLen <= (real_t)0) totalLen = (real_t)1;

    Matrix33r Rshape = Matrix33r::Identity();
    if (section.type == BeamBodyType::Capsule || section.type == BeamBodyType::Cylinder) {
        Rshape = Quaternion4r::FromTwoVectors(Vector3r::UnitZ(), Vector3r::UnitX()).toRotationMatrix();
    }

    entt::entity root = entt::null;
    entt::entity prev = entt::null;
    entt::entity end = entt::null;

    const auto inertial = cardillo::RigidBody::RigidState::inertial();

    bool hasPrevFrame = false;
    Matrix33r prevFrame = Matrix33r::Identity();

    for (const auto& s : samples) {
        const real_t segLen = s.segLen;

        RigidShape shape;
        if (section.type == BeamBodyType::Cube) {
            shape = CubeShape(Vector3r(segLen * (real_t)0.5, section.width * (real_t)0.5, section.height * (real_t)0.5));
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

        const Matrix33r Rlocal = makeStableFrameFromSample(s, hasPrevFrame ? &prevFrame : nullptr);
        prevFrame = Rlocal;
        hasPrevFrame = true;

        Quaternion4r qlocal(Rlocal * Rshape);
        qlocal.normalize();

        cardillo::RigidBody::RigidState segLocal;
        segLocal.position = s.position - splineCOMWorld;
        segLocal.orientation = qlocal;
        segLocal.rotation = qlocal.toRotationMatrix();

        RigidState segState = cardillo::transform::rigidState(segLocal, stateDefaults, inertial);
        segState.position += splineCOMWorld;
        const entt::entity cur = BodyFactory::addRigidBody(sys, shape, segState, segProps);

        cardillo::C_BeamElement be_cur;
        be_cur.l0 = segLen;
        be_cur.l = segLen;
        if (prev != entt::null) {
            be_cur.prev = prev;
            if (!sys.ecs().any_of<cardillo::C_BeamElement>(prev)) {
                cardillo::C_BeamElement be_prev;
                be_prev.l0 = segLen;
                be_prev.l = segLen;
                be_prev.next = cur;
                sys.ecs().emplace<cardillo::C_BeamElement>(prev, be_prev);
            } else {
                sys.ecs().get<cardillo::C_BeamElement>(prev).next = cur;
            }
        }
        sys.ecs().emplace<cardillo::C_BeamElement>(cur, be_cur);

        if (prev != entt::null) {
            ConstraintFactory::addBeamConstraint(sys, prev, cur, springs, section);
            if (collision_mgr) collision_mgr->disablePair(prev, cur);
        }

        if (root == entt::null) root = cur;
        prev = cur;
        end = cur;
    }

    if (loop && root != entt::null && end != entt::null && end != root) {
        ConstraintFactory::addBeamConstraint(sys, end, root, springs, section);
        if (collision_mgr) collision_mgr->disablePair(end, root);
        if (sys.ecs().any_of<cardillo::C_BeamElement>(end)) {
            sys.ecs().get<cardillo::C_BeamElement>(end).next = root;
        }
        if (sys.ecs().any_of<cardillo::C_BeamElement>(root)) {
            sys.ecs().get<cardillo::C_BeamElement>(root).prev = end;
        }
    }

    return {root, end};
}

}  // namespace

std::pair<entt::entity, entt::entity> BeamFactory::createBeam(World& system, const misc::SplinePattern& spline, const BeamCrossSection& section, const BeamSpringParams& springs,
                                                              const RigidState& stateDefaults, const RigidProps& propsDefaults, size_t segments, cardillo::collision::CollisionCoal* collision_mgr) {
    std::vector<BeamSample> samples;
    samples.reserve(segments);

    appendSplineSamples(spline, segments, samples);

    const Vector3r splineCOMWorld = spline.centerOfMass();
    return buildBeamFromSamples(system, samples, spline.isLoop(), section, springs, stateDefaults, propsDefaults, splineCOMWorld, collision_mgr);
}

std::pair<entt::entity, entt::entity> BeamFactory::createBeams(World& system, const std::vector<const misc::SplinePattern*>& splines, const BeamCrossSection& section, const BeamSpringParams& springs,
                                                               const RigidState& stateDefaults, const RigidProps& propsDefaults, size_t segments, cardillo::collision::CollisionCoal* collision_mgr) {
    real_t totalLen = (real_t)0;
    for (const auto* sp : splines) {
        if (sp) totalLen += sp->totalLength();
    }

    if (totalLen <= (real_t)0 || segments == 0) {
        return {entt::null, entt::null};
    }

    std::vector<BeamSample> allSamples;
    allSamples.reserve(segments);

    for (const auto* sp : splines) {
        if (!sp) continue;

        const real_t frac = sp->totalLength() / totalLen;
        const size_t segCount = std::max((size_t)1, (size_t)std::round((real_t)segments * frac));

        appendSplineSamples(*sp, segCount, allSamples);
    }

    if (allSamples.empty()) {
        return {entt::null, entt::null};
    }

    // Treat stitched spline list as one chain so tangent/frame continuity is
    // handled in a single pass without rigid joints between per-spline beams.
    // Preserve previous loop behavior for the single-spline case.
    const bool loop = (splines.size() == 1 && splines[0] != nullptr && splines[0]->isLoop());
    const Vector3r splineCOMWorld = allSamples.front().position;
    return buildBeamFromSamples(system, allSamples, loop, section, springs, stateDefaults, propsDefaults, splineCOMWorld, collision_mgr);
}

}  // namespace cardillo::physics
