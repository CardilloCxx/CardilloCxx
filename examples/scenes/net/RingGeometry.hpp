#pragma once

#include "../SceneBase.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace cardillo;

// Ring geometry shared by the ring net scenes (NetScene, NetDetailedScene).

// Closed ring of centreline radius R in the local x-y plane with the out-of-plane wave
// z = w cos(2 theta), a saddle. Real ring nets are woven: a ring passes over a linked neighbour at
// one crossing and under it at the other. Identical saddle rings whose neighbours lie in the
// diagonal directions (theta = 45°, 135°, ...) do this, so a net of them starts flat and
// interlinked without tilted rings (see saddleAmplitude()).
class SaddleRingSpline : public misc::SplinePattern {
   public:
    SaddleRingSpline(real_t radius, real_t amplitude) : m_r(radius), m_w(amplitude) {
        constexpr int n = 4096;
        for (int i = 0; i < n; ++i) m_len += (sample((real_t)(i + 1) / n).position - sample((real_t)i / n).position).norm();
    }

    real_t totalLength() const override { return m_len; }
    bool isLoop() const override { return true; }
    Vector3r centerOfMass() const override { return Vector3r::Zero(); }

    misc::SplineSample sample(real_t alpha) const override {
        const real_t theta = (real_t)(2 * M_PI) * std::clamp(alpha, (real_t)0, (real_t)1);
        const real_t c = std::cos(theta);
        const real_t s = std::sin(theta);
        const Vector3r position(m_r * c, m_r * s, m_w * std::cos(2 * theta));
        const Vector3r tangent = Vector3r(-m_r * s, m_r * c, -2 * m_w * std::sin(2 * theta)).normalized();
        const Vector3r radial(c, s, 0);  // outward, as the normal of CircleSpline
        const Vector3r normal = (radial - radial.dot(tangent) * tangent).normalized();
        return {position, tangent, normal, tangent.cross(normal)};
    }

   private:
    real_t m_r{(real_t)1};
    real_t m_w{(real_t)0};
    real_t m_len{(real_t)0};
};

// Saddle amplitude w of rings with centreline radius R whose linked neighbour sits at `offset`
// (x-y plane): at the two crossings (angles ±alpha from the line of centres, direction psi) the
// centrelines are 2 w sin(2 alpha) |sin(2 psi)| apart, which is set to `separation`.
inline real_t saddleAmplitude(real_t ringRadius, const Vector3r& offset, real_t separation) {
    const real_t crossingAngle = std::acos(offset.head<2>().norm() / (2 * ringRadius));
    const real_t linkAngle = std::atan2(offset.y(), offset.x());
    return separation / (2 * std::sin(2 * crossingAngle) * std::abs(std::sin(2 * linkAngle)));
}

// Smallest distance between the polygons through the beam nodes (alpha = k / segments) of two
// closed curves shifted by offsetA and offsetB.
inline real_t polygonDistance(const misc::SplinePattern& a, const Vector3r& offsetA, const misc::SplinePattern& b, const Vector3r& offsetB, size_t segments, int subdivisions = 64) {
    auto polygon = [&](const misc::SplinePattern& curve, const Vector3r& offset) {
        std::vector<Vector3r> points;
        for (size_t k = 0; k < segments; ++k) {
            const Vector3r p0 = curve.sample((real_t)k / segments).position + offset;
            const Vector3r p1 = curve.sample((real_t)(k + 1) / segments).position + offset;
            for (int s = 0; s < subdivisions; ++s) points.push_back(p0 + (p1 - p0) * ((real_t)s / subdivisions));
        }
        return points;
    };
    const std::vector<Vector3r> pa = polygon(a, offsetA);
    const std::vector<Vector3r> pb = polygon(b, offsetB);
    real_t d = std::numeric_limits<real_t>::infinity();
    for (const Vector3r& x : pa)
        for (const Vector3r& y : pb) d = std::min(d, (y - x).norm());
    return d;
}
