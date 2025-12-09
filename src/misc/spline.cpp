#include "spline.hpp"
#include <algorithm>

namespace cardillo { namespace misc {

static inline Matrix33r makeFrameFromTangent(const Vector3r &tangent) {
    Vector3r T = tangent.normalized();
    Vector3r up = std::abs(T.z()) < (real_t)0.9 ? Vector3r(0,0,1) : Vector3r(0,1,0);
    Vector3r N = up.cross(T).normalized();
    if (N.squaredNorm() < (real_t)1e-12) {
        up = Vector3r(1,0,0);
        N = up.cross(T).normalized();
    }
    Vector3r B = T.cross(N).normalized();
    Matrix33r M; M.col(0)=T; M.col(1)=N; M.col(2)=B; return M;
}

// LinearSpline ---------------------------------------------------------------
LinearSpline::LinearSpline(const Vector3r &p0, const Vector3r &p1)
: m_p0(p0), m_p1(p1) {
    m_len = (m_p1 - m_p0).norm();
    m_T = (m_len > (real_t)0) ? (m_p1 - m_p0) / m_len : Vector3r(1,0,0);
}

real_t LinearSpline::totalLength() const { return m_len; }

bool LinearSpline::isLoop() const { return false; }

SplineSample LinearSpline::sample(real_t alpha) const {
    alpha = std::min((real_t)1, std::max((real_t)0, alpha));
    Vector3r pos = (1 - alpha) * m_p0 + alpha * m_p1;
    Matrix33r F = makeFrameFromTangent(m_T);
    return { pos, F.col(0), F.col(1), F.col(2) };
}

Vector3r LinearSpline::centerOfMass() const {
    return (m_p0 + m_p1) * (real_t)0.5;
}

// CircleSpline ---------------------------------------------------------------
CircleSpline::CircleSpline(const Vector3r &center, real_t radius, const Vector3r &normal, const Vector3r &dir0, real_t thetaStart, real_t thetaSpan)
: m_c(center), m_r(radius), m_n(normal.normalized()), m_theta0(thetaStart), m_thetaSpan(thetaSpan) {
    Vector3r d = dir0 - dir0.dot(m_n) * m_n; // project onto plane
    if (d.squaredNorm() < (real_t)1e-10) d = Vector3r::UnitX() - Vector3r::UnitX().dot(m_n)*m_n;
    m_u = d.normalized();
    m_v = m_n.cross(m_u).normalized();
}

real_t CircleSpline::totalLength() const { return std::abs(m_thetaSpan) * m_r; }

bool CircleSpline::isLoop() const { return std::abs(std::abs(m_thetaSpan) - (real_t)(2*M_PI)) < (real_t)1e-8; }

SplineSample CircleSpline::sample(real_t alpha) const {
    alpha = std::min((real_t)1, std::max((real_t)0, alpha));
    real_t theta = m_theta0 + m_thetaSpan * alpha;
    Vector3r radial = std::cos(theta) * m_u + std::sin(theta) * m_v; // outward in plane
    Vector3r tangent = (-std::sin(theta) * m_u + std::cos(theta) * m_v).normalized();
    Vector3r normal = radial.normalized();
    Vector3r binormal = tangent.cross(normal).normalized(); // aligns with m_n up to sign
    Vector3r pos = m_c + m_r * radial;
    return { pos, tangent, normal, binormal };
}

Vector3r CircleSpline::centerOfMass() const {
    return m_c;
}

// HelixSpline ----------------------------------------------------------------
HelixSpline::HelixSpline(const Vector3r &center, const Vector3r &axisDir, real_t radius, real_t pitch, real_t turns, const Vector3r &dir0)
: m_c(center), m_d(axisDir.normalized()), m_r(radius), m_pitch(pitch), m_turns(turns) {
    // Build a stable orthonormal basis {u,v} orthogonal to axis m_d using a robust fallback
    Vector3r d = dir0 - dir0.dot(m_d) * m_d; // remove axis component from dir0
    if (d.squaredNorm() < (real_t)1e-12) {
        // dir0 nearly parallel to axis; pick a canonical axis least aligned with m_d
        Vector3r cand = (std::abs(m_d.x()) < (real_t)0.9) ? Vector3r::UnitX() : Vector3r::UnitY();
        d = cand - cand.dot(m_d) * m_d;
    }
    m_u = d.normalized();
    m_v = m_d.cross(m_u).normalized();
}

real_t HelixSpline::totalLength() const {
    real_t perRev = std::sqrt((real_t)(2*M_PI*m_r)*(2*M_PI*m_r) + m_pitch*m_pitch);
    return m_turns * perRev;
}

bool HelixSpline::isLoop() const { return false; }

SplineSample HelixSpline::sample(real_t alpha) const {
    alpha = std::min((real_t)1, std::max((real_t)0, alpha));
    real_t theta = (real_t)2 * M_PI * m_turns * alpha;
    Vector3r radial = std::cos(theta) * m_u + std::sin(theta) * m_v; // outward, |radial|=1
    Vector3r dpos_dtheta = -m_r * std::sin(theta) * m_u + m_r * std::cos(theta) * m_v + (m_pitch / ((real_t)2*M_PI)) * m_d;
    Vector3r tangent = dpos_dtheta.normalized();
    Vector3r normal = radial;
    Vector3r binormal = tangent.cross(normal).normalized();
    Vector3r pos = m_c + m_r * radial + (m_pitch * theta / ((real_t)2*M_PI)) * m_d;
    return { pos, tangent, normal, binormal };
}

Vector3r HelixSpline::centerOfMass() const {
    const real_t Theta = (real_t)2 * M_PI * m_turns;
    if (Theta == (real_t)0) return m_c;
    const real_t sT = std::sin(Theta);
    const real_t cT = std::cos(Theta);
    Vector3r avg_radial = m_r * ((sT / Theta) * m_u + ((1 - cT) / Theta) * m_v);
    Vector3r avg_axial  = ((m_pitch * m_turns) * (real_t)0.5) * m_d; // (pitch*turns)/2
    return m_c + avg_radial + avg_axial;
}

}} // namespace cardillo::misc
