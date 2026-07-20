#pragma once

#include "../SceneBase.hpp"
#include <Eigen/Geometry>

using namespace cardillo;

// Regression test for CondensedAssembler's bilateral elimination-order cache
// (buildBilateralFactorization()): the cache is keyed on the bilateral graph's *structure* (which
// blocks exist, which pairs are coupled), not just observed to be stable in every example scene --
// this scene actually changes that structure at runtime (a spring is added mid-simulation via
// updateScene(), exercising World::markStructureDirty() and DynamicsAssembler::rebuildInteractionW_()'s
// unconditional-every-step rebuild of constraintResults()), so a stale-cache bug here would be
// directly observable rather than just argued from the code.
//
// Layout: three point masses A-B-C in a vertical line. A-B are connected by a spring from the
// start. B-C start unconnected (C free-falls) until t >= addAt, when a second spring connecting
// B-C is added -- growing the bilateral graph from 1 block to 2, changing both `dims` and
// `edgeNodes` as seen by CondensedAssembler::buildBilateralFactorization()'s cache check.
class DynamicConstraintScene : public SceneBase {
public:
    const char* sceneName() const override { return "dynamic_constraint"; }

    void populate(physics::PhysicsEngine& engine) override {
        using namespace cardillo;

        const real_t mass = (real_t)0.2;
        const real_t radius = (real_t)0.03;
        const real_t k = (real_t)500.0;
        const real_t d = (real_t)0.5;

        m_a = engine.addPointMass(mass, Vector3r(0.0, 0.0, 1.0), Vector3r::Zero(), radius);
        m_b = engine.addPointMass(mass, Vector3r(0.0, 0.0, 0.7), Vector3r::Zero(), radius);
        m_c = engine.addPointMass(mass, Vector3r(0.0, 0.0, 0.4), Vector3r::Zero(), radius);

        engine.addLinearDistanceConstraint(m_a, m_b, Vector3r::Zero(), Vector3r::Zero(), k, d);
        m_springK = k;
        m_springD = d;
    }

    void updateScene(physics::PhysicsEngine& engine, real_t t, real_t /*dt*/) override {
        if (!m_added && t >= m_addAt) {
            engine.addLinearDistanceConstraint(m_b, m_c, Vector3r::Zero(), Vector3r::Zero(), m_springK, m_springD);
            m_added = true;
        }
    }

private:
    entt::entity m_a{entt::null}, m_b{entt::null}, m_c{entt::null};
    real_t m_springK{0}, m_springD{0};
    real_t m_addAt{0.02};  // seconds
    bool m_added{false};
};
