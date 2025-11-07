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

// CircleSpline ---------------------------------------------------------------
CircleSpline::CircleSpline(const Vector3r &center, real_t radius, const Vector3r &normal, const Vector3r &dir0)
: m_c(center), m_r(radius), m_n(normal.normalized()) {
    Vector3r d = dir0 - dir0.dot(m_n) * m_n; // project onto plane
    if (d.squaredNorm() < (real_t)1e-10) d = Vector3r::UnitX() - Vector3r::UnitX().dot(m_n)*m_n;
    m_u = d.normalized();
    m_v = m_n.cross(m_u).normalized();
}

real_t CircleSpline::totalLength() const { return (real_t)2 * M_PI * m_r; }

bool CircleSpline::isLoop() const { return true; }

SplineSample CircleSpline::sample(real_t alpha) const {
    alpha = std::min((real_t)1, std::max((real_t)0, alpha));
    real_t theta = (real_t)2 * M_PI * alpha;
    Vector3r radial = std::cos(theta) * m_u + std::sin(theta) * m_v; // outward in plane
    Vector3r tangent = (-std::sin(theta) * m_u + std::cos(theta) * m_v).normalized();
    Vector3r normal = radial.normalized();
    Vector3r binormal = tangent.cross(normal).normalized(); // aligns with m_n up to sign
    Vector3r pos = m_c + m_r * radial;
    return { pos, tangent, normal, binormal };
}

// HelixSpline ----------------------------------------------------------------
HelixSpline::HelixSpline(const Vector3r &center, const Vector3r &axisDir, real_t radius, real_t pitch, real_t turns, const Vector3r &dir0)
: m_c(center), m_d(axisDir.normalized()), m_r(radius), m_pitch(pitch), m_turns(turns) {
    Vector3r d = dir0 - dir0.dot(m_d) * m_d;
    if (d.squaredNorm() < (real_t)1e-10) d = Vector3r::UnitX() - Vector3r::UnitX().dot(m_d)*m_d;
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
    Vector3r radial = std::cos(theta) * m_u + std::sin(theta) * m_v; // outward
    Vector3r dpos_dtheta = -m_r * std::sin(theta) * m_u + m_r * std::cos(theta) * m_v + (m_pitch / ((real_t)2*M_PI)) * m_d;
    Vector3r tangent = dpos_dtheta.normalized();
    Vector3r normal = (radial - radial.dot(tangent)*tangent).normalized();
    Vector3r binormal = tangent.cross(normal).normalized();
    Vector3r pos = m_c + m_r * radial + (m_pitch * theta / ((real_t)2*M_PI)) * m_d;
    return { pos, tangent, normal, binormal };
}

}} // namespace cardillo::misc

namespace cardillo { namespace misc {

std::vector<SplineSample> sampleContinuous(const SplinePattern& spline, size_t count) {
    std::vector<SplineSample> out;
    if (count == 0) return out;
    out.reserve(count);

    // First sample: use provided frame directly
    SplineSample s0 = spline.sample((real_t)0);
    out.push_back(s0);
    Vector3r prevT = s0.tangent;
    Vector3r prevN = s0.normal;
    Vector3r prevB = s0.binormal;

    for (size_t i = 1; i < count; ++i) {
        real_t alpha = (real_t)i / (real_t)(count - 1);
        SplineSample si = spline.sample(alpha);
        Vector3r T = si.tangent.normalized();
        Vector3r axis = prevT.cross(T);
        real_t sinAngle = axis.norm();
        real_t cosAngle = std::max((real_t)-1, std::min((real_t)1, prevT.dot(T)));
        if (cosAngle < (real_t)-0.999 && sinAngle < (real_t)1e-6) {
            prevT = T;
            prevN = -prevN;
            prevB = -prevB;
            prevN = (prevN - prevN.dot(T) * T).normalized();
            prevB = T.cross(prevN).normalized();
        } else if (sinAngle > (real_t)1e-9) {
            axis.normalize();
            real_t angle = std::atan2(sinAngle, cosAngle);
            auto rotateVec = [&](const Vector3r& v){
                return v * std::cos(angle) + axis.cross(v) * std::sin(angle) + axis * (axis.dot(v)) * (1 - std::cos(angle));
            };
            Vector3r Np = rotateVec(prevN);
            Np = (Np - Np.dot(T) * T).normalized();
            Vector3r Bp = T.cross(Np).normalized();
            if (Np.dot(prevN) < (real_t)0) { Np = -Np; Bp = -Bp; }
            prevT = T; prevN = Np; prevB = Bp;
        } else {
            prevT = T;
            prevN = (prevN - prevN.dot(T) * T).normalized();
            prevB = T.cross(prevN).normalized();
        }
        out.push_back({si.position, prevT, prevN, prevB});
    }
    return out;
}

}} // namespace cardillo::misc
