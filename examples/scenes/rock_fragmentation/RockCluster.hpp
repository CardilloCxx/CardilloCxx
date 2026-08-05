#pragma once

#include "RockFragmentationBase.hpp"

// Rock as a cluster of 13 equal spheres on FCC sites (one centre sphere plus its 12 nearest
// neighbours), every pair of touching spheres joined by an elastic bond across a circular section of
// the sphere radius. 36 bonds, so the rock can fall apart into up to 13 fragments. See
// RockFragmentationBase for the bond load measures and the sign convention.
class RockClusterScene : public RockFragmentationBase {
public:
    const char* sceneName() const override { return "rock_cluster"; }

    RockClusterScene() {
        // Calibrated against the peak loads of an unbreakable cluster in this drop. The outcome is
        // sensitive to tau_c: failure of the compressed core bonds is what triggers the cascade.
        m_sigmaC = (real_t)250;
        m_tauC = (real_t)750;
    }

protected:
    real_t supportDistance(const Vector3r& n) const override {
        real_t support = (real_t)0;
        for (const auto& c : sites()) support = std::max(support, -n.dot(c));
        return support + m_radius;
    }

    void buildRock(physics::PhysicsEngine& engine, const Vector3r& center, const Vector3r& v0) override {
        const auto s = sites();
        const physics::RigidProps props = rockProps();

        // Same velocity for every sphere, so the bonds start unloaded
        for (const auto& c : s) {
            m_bodies.push_back(engine.addRigidBody(physics::SphereShape(m_radius), physics::RigidState(center + c, v0), props));
        }

        const real_t A = (real_t)M_PI * m_radius * m_radius;
        const real_t L = (real_t)2 * m_radius;
        const real_t touching = (real_t)2 * m_radius * (real_t)1.001;
        for (size_t i = 0; i < s.size(); ++i) {
            for (size_t j = i + 1; j < s.size(); ++j) {
                const Vector3r d = s[j] - s[i];
                if (d.norm() > touching) continue;
                addBond(engine, i, j, center + (real_t)0.5 * (s[i] + s[j]), d.normalized(), A, m_radius, L);
            }
        }

        std::cout << "[" << m_tag << "] sphere radius " << m_radius << " m, cluster diameter " << 6 * m_radius << " m" << std::endl;
    }

private:
    // FCC sites: centre sphere plus its 12 nearest neighbours at distance 2 r
    std::vector<Vector3r> sites() const {
        std::vector<Vector3r> s;
        s.push_back(Vector3r::Zero());
        const real_t h = std::sqrt((real_t)2) * m_radius;
        for (int axis = 0; axis < 3; ++axis) {
            for (int s1 = -1; s1 <= 1; s1 += 2) {
                for (int s2 = -1; s2 <= 1; s2 += 2) {
                    Vector3r c = Vector3r::Zero();
                    c[axis] = (real_t)s1 * h;
                    c[(axis + 1) % 3] = (real_t)s2 * h;
                    s.push_back(c);
                }
            }
        }
        return s;
    }

    real_t m_radius{(real_t)0.05};  // m, sphere radius; cluster diameter is 6 r
};
