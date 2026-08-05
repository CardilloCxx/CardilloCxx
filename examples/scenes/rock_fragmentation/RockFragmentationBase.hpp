#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>
#include <cmath>
#include <iostream>
#include <vector>

using namespace cardillo;

// Shared frame for the rock fragmentation scenes: a rock is dropped onto an inclined frictional
// plane and falls apart along its internal bonds. Derived scenes differ only in how the rock is
// built -- which bodies it consists of and where the bonds sit:
//   RockCluster  (rock_cluster)  -- 13 sphere primitives on FCC sites, 36 bonds
//   RockQuarters (rock_quarters) -- sphere cut into 4 wedges (meshes), 4 bonds
//
// Bond load per bond, in the bond frame (column 0 of the joint frame is the bond axis, so row 0 is
// normal, rows 1-2 shear, row 3 torsion, rows 4-5 bending):
//   sigma = Lambda_0 / A + |(Lambda_4, Lambda_5)| * R / I      tensile, compression gives sigma < 0
//   tau   = |(Lambda_1, Lambda_2)| / A + |Lambda_3| * R / Jp   shear
// with the bond's section area A and an equivalent circular section radius R (I = pi R^4 / 4,
// Jp = pi R^4 / 2), both supplied per bond by the derived scene. This is the structure of the
// parallel-bond model of Potyondy & Cundall (2004), Int. J. Rock Mech. Min. Sci. 41(8), evaluated on
// constraint impulses (unit N s / m^2) instead of forces, so that the thresholds do not scale with
// the time step. Divide a threshold by dt to read it as a stress.
//
// The sign of Lambda_0 was determined by measurement, not assumed: a sphere hanging from a fixed one
// by a bond whose axis points from the fixed to the hanging sphere gives Lambda_0 = +m g dt, so with
// the bond axis running from body i to body j, positive Lambda_0 is tension.
//
// The bonds must have finite stiffness. Both bond graphs contain cycles, so rigid bonds would impose
// more constraint rows than can be independent: the Schur complement becomes singular (SparseLU
// factorization fails) and the bond loads of the statically indeterminate assembly are not unique,
// which would leave the fracture criterion without a well-defined value. Finite stiffness
// regularizes the system and makes the load sharing between bonds unique.
class RockFragmentationBase : public SceneBase {
public:
    ~RockFragmentationBase() override {
        std::cout << "[" << m_tag << "] " << m_numBroken << " of " << m_bonds.size() << " bonds broken, rock split into " << countFragments() << " fragments" << std::endl;
        std::cout << "[" << m_tag << "] peak bond load: sigma = " << m_peakSigma << " Ns/m^2, tau = " << m_peakTau << " Ns/m^2" << std::endl;
    }

    void populate(physics::PhysicsEngine& engine) override {
        m_tag = sceneName();
        const real_t slope = m_slopeDeg * (real_t)M_PI / (real_t)180.0;
        const Vector3r nGround(std::sin(slope), (real_t)0, std::cos(slope));

        // Plane shifted one half-length downhill inside its own surface, so the rock lands at the
        // upper edge; the shift leaves the halfspace collider unchanged
        const Vector3r up(0, 1, 0);
        const Vector3r across = (up - up.dot(nGround) * nGround).normalized();
        const Vector3r downhill = across.cross(nGround);
        engine.addStaticBody(physics::PlaneShape{nGround, up, m_planeHalfWidth, m_planeHalfLength}, physics::RigidState(downhill * m_planeHalfLength));

        // Rock centre m_dropHeight above the inclined surface, measured along its normal
        const Vector3r center = nGround * (supportDistance(nGround) + m_dropHeight);
        const Vector3r v0(0, 0, -m_initialVelocity);
        buildRock(engine, center, v0);

        real_t mass = 0;
        for (auto e : m_bodies) mass += engine.getMass(e)(0, 0);
        const real_t vImpact = std::sqrt(m_initialVelocity * m_initialVelocity + (real_t)2 * (real_t)9.81 * m_dropHeight / std::cos(slope));
        std::cout << "[" << m_tag << "] " << m_bodies.size() << " bodies, " << m_bonds.size() << " bonds, total mass " << mass << " kg" << std::endl;
        std::cout << "[" << m_tag << "] slope " << m_slopeDeg << " deg, mu = " << m_friction << ", release velocity " << m_initialVelocity << " m/s, impact velocity " << vImpact << " m/s" << std::endl;
        std::cout << "[" << m_tag << "] thresholds: sigma_c = " << m_sigmaC << " Ns/m^2, tau_c = " << m_tauC << " Ns/m^2" << std::endl;
    }

    void updateScene(physics::PhysicsEngine& engine, real_t t, real_t /*dt*/) override {
        // All bond loads are read before any bond is removed: removing one invalidates the row
        // mapping of constraintImpulse() for the rest of this step
        for (auto& b : m_bonds) {
            b.sigma = 0;
            b.tau = 0;
            if (b.broken) continue;

            const VectorXr Lambda = engine.constraintImpulse(b.index);
            if (Lambda.size() < 6) continue;

            b.sigma = Lambda[0] / b.A + Vector2r(Lambda[4], Lambda[5]).norm() * b.R / b.I;
            b.tau = Vector2r(Lambda[1], Lambda[2]).norm() / b.A + std::abs(Lambda[3]) * b.R / b.Jp;
            m_peakSigma = std::max(m_peakSigma, b.sigma);
            m_peakTau = std::max(m_peakTau, b.tau);
        }

        for (auto& b : m_bonds) {
            if (b.broken) continue;

            const bool tensile = b.sigma > m_sigmaC;
            if (!tensile && b.tau <= m_tauC) continue;

            engine.removeConstraint(b.index);
            engine.enableCollisionBetween(m_bodies[b.i], m_bodies[b.j]);
            b.broken = true;
            ++m_numBroken;
            if (m_numBroken == 1) onFirstBreak(engine);
            std::cout << "[" << m_tag << "] bond " << b.i << "-" << b.j << " broken at t = " << t << " s in " << (tensile ? "tension" : "shear") << ": sigma = " << b.sigma << ", tau = " << b.tau
                      << " Ns/m^2" << std::endl;
        }
    }

protected:
    struct Bond {
        size_t index;
        size_t i, j;
        real_t A{0};   // section area, m^2
        real_t R{0};   // equivalent circular section radius, m
        real_t I{0};   // pi R^4 / 4
        real_t Jp{0};  // pi R^4 / 2
        bool broken{false};
        real_t sigma{0}, tau{0};
    };

    // Distance from the rock centre to its lowest surface point along -n
    virtual real_t supportDistance(const Vector3r& n) const = 0;

    // Create the bodies (into m_bodies) and their bonds (via addBond)
    virtual void buildRock(physics::PhysicsEngine& engine, const Vector3r& center, const Vector3r& v0) = 0;

    // Called once, when the first bond breaks
    virtual void onFirstBreak(physics::PhysicsEngine& /*engine*/) {}

    // Elastic 6-DOF bond between bodies i and j across a section of area A and equivalent radius R,
    // with its joint frame at `point` and column 0 along `axis` (pointing from i to j). `length` is
    // the bond length used for the stiffness.
    void addBond(physics::PhysicsEngine& engine, size_t i, size_t j, const Vector3r& point, const Vector3r& axis, real_t A, real_t R, real_t length) {
        const real_t I = (real_t)M_PI * std::pow(R, 4) / (real_t)4;
        const real_t Jp = (real_t)M_PI * std::pow(R, 4) / (real_t)2;
        const Vector3r Ktrans = Vector3r::Constant(m_bondModulus * A / length);
        const Vector3r Krot = Vector3r::Constant(m_bondModulus * I / length);

        const physics::JointFrame frame = physics::JointFrame::fromAxis(point, axis);
        const size_t idx = engine.addTranslationRotationConstraint(m_bodies[i], m_bodies[j], frame, Ktrans, Vector3r::Zero(), Krot, Vector3r::Zero());
        engine.disableCollisionBetween(m_bodies[i], m_bodies[j]);
        m_bonds.push_back(Bond{idx, i, j, A, R, I, Jp});
    }

    physics::RigidProps rockProps() const {
        physics::RigidProps props = physics::RigidProps::withDensity(m_density);
        props.friction = m_friction;
        props.restitution_normal = (real_t)0;
        return props;
    }

    // Connected components of the surviving bond graph
    int countFragments() const {
        std::vector<size_t> root(m_bodies.size());
        for (size_t i = 0; i < root.size(); ++i) root[i] = i;
        const auto find = [&root](size_t x) {
            while (root[x] != x) x = root[x];
            return x;
        };
        for (const auto& b : m_bonds) {
            if (b.broken) continue;
            const size_t ri = find(b.i), rj = find(b.j);
            if (ri != rj) root[ri] = rj;
        }
        int n = 0;
        for (size_t i = 0; i < root.size(); ++i) {
            if (find(i) == i) ++n;
        }
        return n;
    }

    // Material and load case (model hypothesis: granite, rho = 2600 kg/m^3, E = 50 GPa)
    real_t m_density{(real_t)2600};       // kg/m^3
    real_t m_bondModulus{(real_t)50e9};   // Pa, bond Young's modulus
    real_t m_friction{(real_t)0.5};       // Coulomb friction coefficient
    real_t m_dropHeight{(real_t)0.5};     // m, above the inclined surface along its normal
    real_t m_initialVelocity{(real_t)5};  // m/s, downward at release
    real_t m_slopeDeg{(real_t)30};        // deg, inclination of the plane
    real_t m_planeHalfWidth{(real_t)2};   // m, visual half-width of the plane across the slope
    real_t m_planeHalfLength{(real_t)4};  // m, visual half-length of the plane along the slope

    // Bond strengths, set by the derived scene (calibrated per rock)
    real_t m_sigmaC{(real_t)1e30};  // Ns/m^2, tensile
    real_t m_tauC{(real_t)1e30};    // Ns/m^2, shear

    std::vector<entt::entity> m_bodies;
    std::vector<Bond> m_bonds;
    int m_numBroken{0};
    real_t m_peakSigma{0}, m_peakTau{0};
    std::string m_tag{"rock"};
};
