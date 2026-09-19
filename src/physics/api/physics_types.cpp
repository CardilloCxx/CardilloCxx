#include "physics_types.hpp"

#include "../../rigid_body/transformations.hpp"
#include "../world.hpp"

namespace cardillo::physics {

RigidState::RigidState(const Vector3r& p_local, const Vector3r& v_local, const Quaternion4r& q_local, const Vector3r& w_local, entt::entity refEntity, entt::registry& reg) {
    const RigidBody::RigidState refState = RigidBody::getState(reg, refEntity);
    const RigidBody::RigidState inertial = RigidBody::RigidState::inertial();

    RigidState localState;
    localState.position = p_local;
    localState.setOrientation(q_local);
    localState.linearVelocity = v_local;
    localState.angularVelocity = w_local;

    *this = transform::rigidState(localState, refState, inertial);
}

Vector2r BeamCrossSection::centroidOf(const std::vector<Vector2r>& poly) {
    real_t A = 0, Cx = 0, Cy = 0;
    size_t n = poly.size();
    for (size_t i = 0; i < n; ++i) {
        const Vector2r& p1 = poly[i];
        const Vector2r& p2 = poly[(i + 1) % n];
        real_t cross = p1.x() * p2.y() - p2.x() * p1.y();
        A += cross;
        Cx += (p1.x() + p2.x()) * cross;
        Cy += (p1.y() + p2.y()) * cross;
    }
    A *= (real_t)0.5;
    if (std::abs(A) < (real_t)1e-12) {
        Vector2r avg(0, 0);
        for (const auto& p : poly) avg += p;
        return n > 0 ? avg / (real_t)n : avg;
    }
    return Vector2r(Cx / ((real_t)6.0 * A), Cy / ((real_t)6.0 * A));
}

std::vector<Vector2r> BeamCrossSection::recenter(const std::vector<Vector2r>& poly) {
    Vector2r c = centroidOf(poly);
    std::vector<Vector2r> out;
    out.reserve(poly.size());
    for (const auto& p : poly) out.push_back(p - c);
    return out;
}

BeamCrossSection BeamCrossSection::square(real_t w, real_t h) {
    BeamCrossSection sec;
    sec.width = w;
    sec.height = h;
    sec.type = BeamCrossSectionType::Square;
    sec.polygon = recenter({Vector2r(-w / 2, -h / 2), Vector2r(w / 2, -h / 2),
                            Vector2r(w / 2, h / 2), Vector2r(-w / 2, h / 2)});
    return sec;
}

BeamCrossSection BeamCrossSection::triangle(real_t w, real_t h) {
    BeamCrossSection sec;
    sec.width = w;
    sec.height = h;
    sec.type = BeamCrossSectionType::Triangle;
    sec.polygon = recenter({Vector2r(-w / 2, -h / 2), Vector2r(w / 2, -h / 2), Vector2r(0, h / 2)});
    return sec;
}

BeamCrossSection BeamCrossSection::round(real_t radius, size_t numSegments) {
    BeamCrossSection sec;
    sec.radius = radius;
    sec.type = BeamCrossSectionType::Round;
    sec.polygon.reserve(numSegments);
    for (size_t i = 0; i < numSegments; ++i) {
        real_t angle = (real_t)i / (real_t)numSegments * (real_t)2.0 * (real_t)M_PI;
        sec.polygon.push_back(Vector2r(radius * std::cos(angle), radius * std::sin(angle)));
    }
    sec.polygon = recenter(sec.polygon);
    return sec;
}

BeamCrossSection BeamCrossSection::custom(const std::vector<Vector2r>& poly) {
    BeamCrossSection sec;
    sec.polygon = recenter(poly);
    sec.type = BeamCrossSectionType::Polygon;
    return sec;
}

real_t BeamCrossSection::area() const {
    switch (type) {
        case BeamCrossSectionType::Round:
            return (real_t)M_PI * radius * radius;
        case BeamCrossSectionType::Triangle:
            return (real_t)0.5 * width * height;
        case BeamCrossSectionType::Polygon: {
            real_t a = (real_t)0;
            size_t n = polygon.size();
            for (size_t i = 0; i < n; ++i) {
                const Vector2r& p1 = polygon[i];
                const Vector2r& p2 = polygon[(i + 1) % n];
                a += p1.x() * p2.y() - p2.x() * p1.y();
            }
            return std::abs(a) * (real_t)0.5;
        }
        default:
            return width * height;
    }
}

real_t BeamCrossSection::Iy() const {
    switch (type) {
        case BeamCrossSectionType::Round:
            return (real_t)M_PI * std::pow(radius, 4) / (real_t)4.0;
        case BeamCrossSectionType::Triangle:
            return width * std::pow(height, (real_t)3) / (real_t)36.0;
        case BeamCrossSectionType::Polygon: {
            real_t Iy_val = (real_t)0;
            size_t n = polygon.size();
            for (size_t i = 0; i < n; ++i) {
                const Vector2r& p1 = polygon[i];
                const Vector2r& p2 = polygon[(i + 1) % n];
                Iy_val += (p1.x() * p2.y() - p2.x() * p1.y()) * (p1.y() * p1.y() + p1.y() * p2.y() + p2.y() * p2.y());
            }
            return std::abs(Iy_val) / (real_t)12.0;
        }
        default:
            return width * std::pow(height, (real_t)3) / (real_t)12.0;
    }
}

real_t BeamCrossSection::Iz() const {
    switch (type) {
        case BeamCrossSectionType::Round:
            return (real_t)M_PI * std::pow(radius, 4) / (real_t)4.0;
        case BeamCrossSectionType::Triangle:
            return std::pow(width, (real_t)3) * height / (real_t)36.0;
        case BeamCrossSectionType::Polygon: {
            real_t Iz_val = (real_t)0;
            size_t n = polygon.size();
            for (size_t i = 0; i < n; ++i) {
                const Vector2r& p1 = polygon[i];
                const Vector2r& p2 = polygon[(i + 1) % n];
                Iz_val += (p1.x() * p2.y() - p2.x() * p1.y()) * (p1.x() * p1.x() + p1.x() * p2.x() + p2.x() * p2.x());
            }
            return std::abs(Iz_val) / (real_t)12.0;
        }
        default:
            return std::pow(width, (real_t)3) * height / (real_t)12.0;
    }
}

real_t BeamCrossSection::Jp() const { return Iy() + Iz(); }

real_t BeamCrossSection::maxAbsY() const {
    real_t m = (real_t)0;
    for (const auto& p : polygon) m = std::max(m, std::abs(p.y()));
    return m;
}

real_t BeamCrossSection::maxAbsX() const {
    real_t m = (real_t)0;
    for (const auto& p : polygon) m = std::max(m, std::abs(p.x()));
    return m;
}

real_t BeamCrossSection::sectionModulus() const {
    switch (type) {
        case BeamCrossSectionType::Round:
            return (real_t)M_PI * std::pow(radius, 3) / (real_t)4.0;
        case BeamCrossSectionType::Triangle: {
            real_t S_horizontal = width * height * height / (real_t)24.0;
            real_t S_vertical   = width * width * height / (real_t)18.0;
            return std::min(S_horizontal, S_vertical);
        }
        case BeamCrossSectionType::Polygon: {
            real_t cy = maxAbsY();
            real_t cz = maxAbsX();
            real_t Wy = (cy > (real_t)0) ? Iy() / cy : (real_t)0;
            real_t Wz = (cz > (real_t)0) ? Iz() / cz : (real_t)0;
            return std::min(Wy, Wz);
        }
        default:
            return std::max(width, height) * std::min(width, height) * std::min(width, height) / (real_t)6.0;
    }
}

BeamHullShape::BeamHullShape(const BeamCrossSection& cs, float len) : cross_section(cs), length(len) {}

BeamSpringParams::BeamSpringParams(const Vector3r& Ke_in, const Vector3r& Kf_in, real_t dampingFactor_in)
    : Ke_direct(Ke_in), Kf_direct(Kf_in), dampingFactor(dampingFactor_in) {}

Vector3r BeamSpringParams::Ke(real_t segLen, const BeamCrossSection& sec) const {
    if (Ke_direct.has_value()) return *Ke_direct;
    const real_t G = E / ((real_t)2.0 * ((real_t)1.0 + nu));
    const real_t A = sec.area();
    Vector3r base(E * A / segLen, G * A / segLen, G * A / segLen);
    return base.cwiseProduct(scaleKe);
}

Vector3r BeamSpringParams::Kf(real_t segLen, const BeamCrossSection& sec) const {
    if (Kf_direct.has_value()) return *Kf_direct;
    const real_t G = E / ((real_t)2.0 * ((real_t)1.0 + nu));
    Vector3r base(G * sec.Jp() / segLen, E * sec.Iy() / segLen, E * sec.Iz() / segLen);
    return base.cwiseProduct(scaleKf);
}

void BeamSpringParams::setDampingFromFactor(real_t d) { dampingFactor = d; }

BeamSpringParams BeamSpringParams::fromMaterial(real_t E_in, real_t nu_in, real_t axialScale, real_t shearScale, 
                                                 real_t torsionScale, real_t bendYScale, real_t bendZScale, 
                                                 real_t dampingFactor_in) {
    BeamSpringParams p;
    p.E = E_in;
    p.nu = nu_in;
    p.scaleKe = Vector3r(axialScale, shearScale, shearScale);
    p.scaleKf = Vector3r(torsionScale, bendYScale, bendZScale);
    p.dampingFactor = dampingFactor_in;
    return p;
}

}  // namespace cardillo::physics
