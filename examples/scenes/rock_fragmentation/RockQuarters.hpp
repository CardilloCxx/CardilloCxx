#pragma once

#include "RockFragmentationBase.hpp"

// Rock as a sphere cut by two orthogonal planes through the polar axis into four equal wedges, each
// loaded as a mesh (res/meshes/sphere_quarter.obj). Neighbouring wedges are joined by an elastic bond
// across their shared flat interface, a half-disc of the sphere radius: A = pi R^2 / 2 exactly, while
// the lever arms use a circular section of equal area (approximation), R_eq = R / sqrt(2). The exact
// second moments of a half-disc are direction dependent (about 0.110 R^4 and 0.393 R^4 about its two
// centroidal axes), the equal-area circle (0.196 R^4) lies between them. See RockFragmentationBase
// for the bond load measures and the sign convention.
//
// The four bonds form a cycle, so the rock splits into 2 fragments if two opposite bonds fail and
// into 4 if all of them do.
class RockQuartersScene : public RockFragmentationBase {
public:
    const char* sceneName() const override { return "rock_quarters"; }

    RockQuartersScene() {
        // Calibrated against the peak loads of an unbreakable rock in this drop (sigma = 1330,
        // tau = 998 Ns/m^2)
        m_sigmaC = (real_t)500;
        m_tauC = (real_t)600;
    }

protected:
    real_t supportDistance(const Vector3r& /*n*/) const override { return m_radius; }

    void buildRock(physics::PhysicsEngine& engine, const Vector3r& center, const Vector3r& v0) override {
        const physics::RigidProps props = rockProps();

        // The mesh origin is the sphere centre; the four wedges differ only by a rotation about z
        const physics::MeshShape wedge(std::string(PROJECT_SOURCE_DIR) + "/res/meshes/sphere_quarter.obj", Vector3r::Constant(m_radius));
        for (int i = 0; i < 4; ++i) {
            const Quaternion4r q(Eigen::AngleAxis<real_t>((real_t)i * (real_t)M_PI_2, Vector3r::UnitZ()));
            m_bodies.push_back(engine.addRigidBody(wedge, physics::RigidState(center, v0, q), props));
        }

        // Interface between wedge i and i+1: half-disc in the plane at the middle of the gap between
        // them, centroid at 4 R / (3 pi) from the sphere centre, normal along +e_phi
        const real_t A = (real_t)0.5 * (real_t)M_PI * m_radius * m_radius;
        const real_t Req = m_radius / std::sqrt((real_t)2);
        const real_t centroid = (real_t)4 * m_radius / ((real_t)3 * (real_t)M_PI);

        for (size_t i = 0; i < m_bodies.size(); ++i) {
            const size_t j = (i + 1) % m_bodies.size();
            const real_t phi = (real_t)(i + 1) * (real_t)M_PI_2 - (real_t)0.5 * m_gapAngle * (real_t)M_PI / (real_t)180.0;
            const Vector3r er(std::cos(phi), std::sin(phi), (real_t)0);
            const Vector3r ephi(-std::sin(phi), std::cos(phi), (real_t)0);

            // Bond length: distance between the two wedge centres of mass
            const real_t L = (engine.getPosition(m_bodies[j]).head<3>() - engine.getPosition(m_bodies[i]).head<3>()).norm();
            addBond(engine, i, j, center + centroid * er, ephi, A, Req, L);

            if (i == 0) std::cout << "[" << m_tag << "] sphere diameter " << 2 * m_radius << " m, bond length " << L << " m, interface area " << A << " m^2" << std::endl;
        }

        // The two opposite pairs share only the polar axis, where all four wedges meet: keep them
        // apart as long as the rock is intact (addBond already disabled the four bonded pairs)
        engine.disableCollisionBetween(m_bodies[0], m_bodies[2]);
        engine.disableCollisionBetween(m_bodies[1], m_bodies[3]);
    }

    void onFirstBreak(physics::PhysicsEngine& engine) override {
        engine.enableCollisionBetween(m_bodies[0], m_bodies[2]);
        engine.enableCollisionBetween(m_bodies[1], m_bodies[3]);
    }

private:
    real_t m_radius{(real_t)0.15};   // m, sphere radius
    real_t m_gapAngle{(real_t)0.4};  // deg, azimuthal gap between neighbouring wedges
};
