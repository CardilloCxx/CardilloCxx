#include "spline.hpp"
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <limits>

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

// CatmullRomSpline -----------------------------------------------------------
CatmullRomSpline::CatmullRomSpline(std::vector<Vector3r> controlPoints, bool loop)
: m_cp(std::move(controlPoints)), m_loop(loop) {
    if (m_cp.size() < 2) {
        m_len = (real_t)0;
    } else {
        buildArcLengthLUT(20);
    }
}

real_t CatmullRomSpline::totalLength() const { return m_len; }

bool CatmullRomSpline::isLoop() const { return m_loop; }

Vector3r CatmullRomSpline::controlPoint(int idx) const {
    const int n = static_cast<int>(m_cp.size());
    if (n == 0) return Vector3r::Zero();
    if (m_loop) {
        int j = idx % n;
        if (j < 0) j += n;
        return m_cp[j];
    }
    if (idx < 0) return m_cp.front();
    if (idx >= n) return m_cp.back();
    return m_cp[static_cast<size_t>(idx)];
}

Vector3r CatmullRomSpline::segmentPosition(int seg, real_t t) const {
    const Vector3r p0 = controlPoint(seg - 1);
    const Vector3r p1 = controlPoint(seg);
    const Vector3r p2 = controlPoint(seg + 1);
    const Vector3r p3 = controlPoint(seg + 2);

    const real_t t2 = t * t;
    const real_t t3 = t2 * t;
    return (real_t)0.5 * ((real_t)2 * p1
        + (-p0 + p2) * t
        + ((real_t)2 * p0 - (real_t)5 * p1 + (real_t)4 * p2 - p3) * t2
        + (-p0 + (real_t)3 * p1 - (real_t)3 * p2 + p3) * t3);
}

Vector3r CatmullRomSpline::segmentTangent(int seg, real_t t) const {
    const Vector3r p0 = controlPoint(seg - 1);
    const Vector3r p1 = controlPoint(seg);
    const Vector3r p2 = controlPoint(seg + 1);
    const Vector3r p3 = controlPoint(seg + 2);

    const real_t t2 = t * t;
    Vector3r deriv = (real_t)0.5 * ((-p0 + p2)
        + ((real_t)2 * ((real_t)2 * p0 - (real_t)5 * p1 + (real_t)4 * p2 - p3)) * t
        + ((real_t)3 * (-p0 + (real_t)3 * p1 - (real_t)3 * p2 + p3)) * t2);
    const real_t nrm = deriv.norm();
    if (nrm <= (real_t)0) return Vector3r::UnitX();
    return deriv / nrm;
}

SplineSample CatmullRomSpline::sample(real_t alpha) const {
    if (m_cp.empty()) {
        Matrix33r F = makeFrameFromTangent(Vector3r::UnitX());
        return { Vector3r::Zero(), F.col(0), F.col(1), F.col(2) };
    }
    if (m_cp.size() == 1) {
        Matrix33r F = makeFrameFromTangent(Vector3r::UnitX());
        return { m_cp.front(), F.col(0), F.col(1), F.col(2) };
    }

    alpha = std::min((real_t)1, std::max((real_t)0, alpha));
    const int segCount = m_loop ? static_cast<int>(m_cp.size()) : static_cast<int>(m_cp.size()) - 1;
    if (segCount <= 0) {
        Matrix33r F = makeFrameFromTangent(Vector3r::UnitX());
        return { m_cp.front(), F.col(0), F.col(1), F.col(2) };
    }

    real_t u = alpha * (real_t)segCount;
    if (!m_arcLut.empty() && m_len > (real_t)0 && m_lutSamplesPerSeg > 0 && m_lutSegCount == segCount) {
        real_t target = alpha * m_len;
        auto it = std::lower_bound(m_arcLut.begin(), m_arcLut.end(), target);
        if (it == m_arcLut.begin()) {
            u = (real_t)0;
        } else if (it == m_arcLut.end()) {
            u = (real_t)segCount;
        } else {
            const std::size_t idx = static_cast<std::size_t>(std::distance(m_arcLut.begin(), it));
            const real_t s1 = m_arcLut[idx - 1];
            const real_t s2 = m_arcLut[idx];
            const real_t denom = s2 - s1;
            const real_t f = (denom > (real_t)0) ? (target - s1) / denom : (real_t)0;
            const real_t sampleIndex = (real_t)(idx - 1) + f;
            u = sampleIndex / (real_t)m_lutSamplesPerSeg;
        }
    }
    int seg = static_cast<int>(std::floor(u));
    real_t t = u - (real_t)seg;
    if (seg >= segCount) {
        seg = segCount - 1;
        t = (real_t)1;
    }

    Vector3r pos = segmentPosition(seg, t);
    Vector3r tangent = segmentTangent(seg, t);
    Matrix33r F = makeFrameFromTangent(tangent);
    return { pos, F.col(0), F.col(1), F.col(2) };
}

real_t CatmullRomSpline::approximateLength(int samplesPerSegment) const {
    if (m_cp.size() < 2) return (real_t)0;
    const int segCount = m_loop ? static_cast<int>(m_cp.size()) : static_cast<int>(m_cp.size()) - 1;
    if (segCount <= 0 || samplesPerSegment <= 0) return (real_t)0;

    real_t length = (real_t)0;
    for (int s = 0; s < segCount; ++s) {
        Vector3r prev = segmentPosition(s, (real_t)0);
        for (int i = 1; i <= samplesPerSegment; ++i) {
            real_t t = (real_t)i / (real_t)samplesPerSegment;
            Vector3r cur = segmentPosition(s, t);
            length += (cur - prev).norm();
            prev = cur;
        }
    }
    return length;
}

void CatmullRomSpline::buildArcLengthLUT(int samplesPerSegment) {
    m_arcLut.clear();
    m_lutSamplesPerSeg = 0;
    m_lutSegCount = 0;

    if (m_cp.size() < 2 || samplesPerSegment <= 0) {
        m_len = (real_t)0;
        return;
    }

    const int segCount = m_loop ? static_cast<int>(m_cp.size()) : static_cast<int>(m_cp.size()) - 1;
    if (segCount <= 0) {
        m_len = (real_t)0;
        return;
    }

    const int totalSamples = segCount * samplesPerSegment;
    m_arcLut.resize(static_cast<std::size_t>(totalSamples) + 1);
    m_arcLut[0] = (real_t)0;

    real_t length = (real_t)0;
    Vector3r prev = segmentPosition(0, (real_t)0);
    int idx = 0;
    for (int s = 0; s < segCount; ++s) {
        for (int i = 1; i <= samplesPerSegment; ++i) {
            real_t t = (real_t)i / (real_t)samplesPerSegment;
            Vector3r cur = segmentPosition(s, t);
            length += (cur - prev).norm();
            ++idx;
            m_arcLut[static_cast<std::size_t>(idx)] = length;
            prev = cur;
        }
    }

    m_len = length;
    m_lutSamplesPerSeg = samplesPerSegment;
    m_lutSegCount = segCount;
}

Vector3r CatmullRomSpline::centerOfMass() const {
    if (m_cp.empty()) return Vector3r::Zero();
    if (m_cp.size() == 1) return m_cp.front();
    const int samples = 100;
    Vector3r accum = Vector3r::Zero();
    real_t total = (real_t)0;
    Vector3r prev = sample((real_t)0).position;
    for (int i = 1; i <= samples; ++i) {
        real_t a = (real_t)i / (real_t)samples;
        Vector3r cur = sample(a).position;
        real_t segLen = (cur - prev).norm();
        accum += (prev + cur) * ((real_t)0.5 * segLen);
        total += segLen;
        prev = cur;
    }
    if (total <= (real_t)0) return m_cp.front();
    return accum / total;
}

// BCC loader ----------------------------------------------------------------
#pragma pack(push, 1)
struct BCCHeader {
    char sign[3];
    unsigned char byteCount;
    char curveType[2];
    char dimensions;
    char upDimension;
    std::uint64_t curveCount;
    std::uint64_t totalControlPointCount;
    char fileInfo[40];
};
#pragma pack(pop)

std::vector<std::shared_ptr<SplinePattern>> loadSplinesFromBCC(const std::string &filePath,
                                                               real_t scale) {
    std::vector<std::shared_ptr<SplinePattern>> splines;
    std::ifstream in(filePath, std::ios::binary);
    if (!in) return splines;

    BCCHeader header{};
    in.read(reinterpret_cast<char*>(&header), sizeof(header));
    if (!in) return splines;

    if (header.sign[0] != 'B' || header.sign[1] != 'C' || header.sign[2] != 'C') return splines;
    if (header.byteCount != 0x44) return splines; // Only supporting 4-byte integers and floats
    if (header.curveType[0] != 'C' || header.curveType[1] != '0') return splines; // Catmull-Rom, uniform
    if (header.dimensions != 3) return splines;

    if (header.curveCount > (std::uint64_t)std::numeric_limits<std::size_t>::max()) return splines;
    splines.reserve(static_cast<std::size_t>(header.curveCount));

    for (std::uint64_t i = 0; i < header.curveCount; ++i) {
        std::int32_t curveControlPointCount = 0;
        in.read(reinterpret_cast<char*>(&curveControlPointCount), sizeof(curveControlPointCount));
        if (!in) break;

        bool isLoop = curveControlPointCount < 0;
        if (curveControlPointCount < 0) curveControlPointCount = -curveControlPointCount;
        if (curveControlPointCount < 0) break;

        const std::size_t count = static_cast<std::size_t>(curveControlPointCount);
        std::vector<Vector3r> cps;
        cps.reserve(count);

        std::vector<float> raw;
        raw.resize(count * 3);
        if (!raw.empty()) {
            in.read(reinterpret_cast<char*>(raw.data()), sizeof(float) * raw.size());
            if (!in) break;
            for (std::size_t j = 0; j < count; ++j) {
                cps.emplace_back((real_t)raw[3*j + 0] * scale,
                                 (real_t)raw[3*j + 1] * scale,
                                 (real_t)raw[3*j + 2] * scale);
            }
        }

        splines.emplace_back(std::make_shared<CatmullRomSpline>(std::move(cps), isLoop));
    }

    return splines;
}

}} // namespace cardillo::misc
