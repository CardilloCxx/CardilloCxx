#pragma once

#include <Eigen/Geometry>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include "types.hpp"

namespace cardillo {
namespace misc {

// Spline sample result: geometric + differential frame at a parameter alpha
struct SplineSample {
    Vector3r position;  // world-space position
    Vector3r tangent;   // unit tangent (direction of increasing alpha)
    Vector3r normal;    // principal normal (orthogonal to tangent)
    Vector3r binormal;  // binormal = tangent x normal (right-handed)
};

// Base spline pattern: alpha in [0,1]
class SplinePattern {
   public:
    virtual ~SplinePattern() = default;
    virtual real_t totalLength() const = 0;
    virtual bool isLoop() const = 0;
    // Sample full geometric/differential data at alpha (clamps to [0,1])
    virtual SplineSample sample(real_t alpha) const = 0;
    virtual Vector3r centerOfMass() const = 0;
};

// Concrete splines (definitions in spline.cpp)
class LinearSpline : public SplinePattern {
   public:
    LinearSpline(const Vector3r& p0, const Vector3r& p1);
    real_t totalLength() const override;
    bool isLoop() const override;
    SplineSample sample(real_t alpha) const override;
    Vector3r centerOfMass() const override;

   private:
    Vector3r m_p0{Vector3r::Zero()}, m_p1{Vector3r::Zero()};
    real_t m_len{(real_t)0};
    Vector3r m_T{Vector3r::UnitX()};
};

class CircleSpline : public SplinePattern {
   public:
    // Full circle by default (theta in [0, 2pi)). Optionally specify start angle and span
    // (radians).
    CircleSpline(const Vector3r& center, real_t radius, const Vector3r& normal = Vector3r::UnitZ(), const Vector3r& dir0 = Vector3r::UnitX(), real_t thetaStart = (real_t)0,
                 real_t thetaSpan = (real_t)(2 * M_PI));
    real_t totalLength() const override;
    bool isLoop() const override;
    SplineSample sample(real_t alpha) const override;
    Vector3r centerOfMass() const override;

   private:
    Vector3r m_c{Vector3r::Zero()};
    real_t m_r{(real_t)1};
    Vector3r m_n{Vector3r::UnitZ()};
    Vector3r m_u{Vector3r::UnitX()};
    Vector3r m_v{Vector3r::UnitY()};
    real_t m_theta0{(real_t)0};
    real_t m_thetaSpan{(real_t)(2 * M_PI)};
};

class HelixSpline : public SplinePattern {
   public:
    HelixSpline(const Vector3r& center, const Vector3r& axisDir, real_t radius, real_t pitch, real_t turns, const Vector3r& dir0 = Vector3r::UnitX());
    real_t totalLength() const override;
    bool isLoop() const override;
    SplineSample sample(real_t alpha) const override;
    Vector3r centerOfMass() const override;

   private:
    Vector3r m_c{Vector3r::Zero()};
    Vector3r m_d{Vector3r::UnitZ()};
    real_t m_r{(real_t)1};
    real_t m_pitch{(real_t)0.1};
    real_t m_turns{(real_t)1};
    Vector3r m_u{Vector3r::UnitX()};
    Vector3r m_v{Vector3r::UnitY()};
};

// Catmull-Rom spline (uniform parameterization, supports open or closed loops)
class CatmullRomSpline : public SplinePattern {
   public:
    CatmullRomSpline(std::vector<Vector3r> controlPoints, bool loop);
    real_t totalLength() const override;
    bool isLoop() const override;
    SplineSample sample(real_t alpha) const override;
    Vector3r centerOfMass() const override;

   private:
    Vector3r controlPoint(int idx) const;
    Vector3r segmentPosition(int seg, real_t t) const;
    Vector3r segmentTangent(int seg, real_t t) const;
    real_t approximateLength(int samplesPerSegment) const;
    void buildArcLengthLUT(int samplesPerSegment);

    std::vector<Vector3r> m_cp;
    bool m_loop{false};
    real_t m_len{(real_t)0};
    std::vector<real_t> m_arcLut;
    int m_lutSamplesPerSeg{0};
    int m_lutSegCount{0};
};

// Load a list of splines from a BCC file (Catmull-Rom)
std::vector<std::shared_ptr<SplinePattern>> loadSplinesFromBCC(const std::string& filePath, real_t scale = (real_t)1);

}  // namespace misc
}  // namespace cardillo
