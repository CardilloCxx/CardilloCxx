#include "derived_entity_sync.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "../../rigid_body/rigid_body.hpp"
#include "trajectory.hpp"

namespace cardillo {
namespace physics {

void DerivedEntitySync::updateBeamElementEntity(World& world, entt::entity e) {
    auto& reg = world.ecs();

    if (!reg.valid(e) || !reg.any_of<cardillo::C_BeamElement, cardillo::C_Position3>(e)) return;
    auto& be = reg.get<cardillo::C_BeamElement>(e);
    real_t newLen = be.l0;

    auto getDesiredLengthBetween = [&](entt::entity a, entt::entity b) -> real_t {
        if (!(reg.valid(a) && reg.valid(b))) {
            std::cout << "Beam length: invalid entities a=" << (int)a << " b=" << (int)b << "\n";
        }
        if (!reg.any_of<cardillo::C_Position3>(a) || !reg.any_of<cardillo::C_Position3>(b)) {
            std::cout << "Beam length: missing positions a=" << reg.any_of<cardillo::C_Position3>(a) << " b=" << reg.any_of<cardillo::C_Position3>(b) << " for entities a=" << (int)a << " b=" << (int)b
                      << "\n";
        }
        if (!reg.any_of<cardillo::C_Orientation>(a) || !reg.any_of<cardillo::C_Orientation>(b)) {
            std::cout << "Beam length: missing orientations a=" << reg.any_of<cardillo::C_Orientation>(a) << " b=" << reg.any_of<cardillo::C_Orientation>(b) << " for entities a=" << (int)a
                      << " b=" << (int)b << "\n";
        }
        if (reg.valid(a) && reg.valid(b) && reg.any_of<cardillo::C_Position3>(a) && reg.any_of<cardillo::C_Position3>(b) && reg.any_of<cardillo::C_Orientation>(a) &&
            reg.any_of<cardillo::C_Orientation>(b)) {
            const auto& pa = reg.get<cardillo::C_Position3>(a).value;
            const auto& pb = reg.get<cardillo::C_Position3>(b).value;
            auto r_AB = pb - pa;
            auto R_A = reg.get<cardillo::C_Orientation>(a).value.toRotationMatrix();
            auto R_B = reg.get<cardillo::C_Orientation>(b).value.toRotationMatrix();

            index_t x_col_A = (reg.any_of<cardillo::C_Capsule>(a) || reg.any_of<cardillo::C_Cylinder>(a)) ? 2 : 0;
            index_t x_col_B = (reg.any_of<cardillo::C_Capsule>(b) || reg.any_of<cardillo::C_Cylinder>(b)) ? 2 : 0;

            auto e_Ax = R_A.col(x_col_A);
            auto e_Bx = R_B.col(x_col_B);
            real_t base = r_AB.norm();
            if (base > (real_t)0) {
                auto r_hat = r_AB.normalized();
                real_t d1 = std::clamp(e_Ax.dot(r_hat), (real_t)-1, (real_t)1);
                real_t d2 = std::clamp(e_Bx.dot(r_hat), (real_t)-1, (real_t)1);
                real_t s1 = std::sqrt(std::max((real_t)0, (real_t)1 - d1 * d1));
                real_t s2 = std::sqrt(std::max((real_t)0, (real_t)1 - d2 * d2));
                real_t base_angle = std::atan2(s1 + s2, d1 + d2);
                real_t c = std::cos(base_angle);
                if (std::abs(c) > (real_t)1e-12) return base / c;
            }
        }
        return (real_t)0;
    };

    index_t contributions = 0;
    real_t totalLen = (real_t)0;

    if (be.prev.has_value()) {
        totalLen += getDesiredLengthBetween(be.prev.value(), e);
        contributions++;
    }

    if (be.next.has_value()) {
        totalLen += getDesiredLengthBetween(e, be.next.value());
        contributions++;
    }

    if (contributions > 0) {
        real_t avgLen = totalLen / (real_t)contributions;
        newLen = avgLen;
    }

    const real_t prevLen = be.l;
    be.l = newLen;

    const real_t eps = (real_t)1e-8;
    if (std::abs(be.l - prevLen) > eps) {
        bool shapeChanged = false;
        if (reg.any_of<cardillo::C_Cube>(e)) {
            auto& cb = reg.get<cardillo::C_Cube>(e);
            const real_t newHalfX = be.l * (real_t)0.5;
            if (std::abs(cb.halfExtents.x() - newHalfX) > eps) {
                cb.halfExtents.x() = newHalfX;
                shapeChanged = true;
            }
        }
        if (reg.any_of<cardillo::C_Capsule>(e)) {
            auto& cap = reg.get<cardillo::C_Capsule>(e);
            const real_t newHalf = be.l * (real_t)0.5;
            if (std::abs(cap.halfLength - newHalf) > eps) {
                cap.halfLength = newHalf;
                shapeChanged = true;
            }
        }
        if (reg.any_of<cardillo::C_Cylinder>(e)) {
            auto& cyl = reg.get<cardillo::C_Cylinder>(e);
            const real_t newHalf = be.l * (real_t)0.5;
            if (std::abs(cyl.halfLength - newHalf) > eps) {
                cyl.halfLength = newHalf;
                shapeChanged = true;
            }
        }
        if (reg.any_of<cardillo::C_RB_Cube>(e)) {
            auto& cb = reg.get<cardillo::C_RB_Cube>(e);
            const real_t newHalfX = be.l * (real_t)0.5;
            if (std::abs(cb.halfExtents.x() - newHalfX) > eps) {
                cb.halfExtents.x() = newHalfX;
                shapeChanged = true;
            }
        }
        if (reg.any_of<cardillo::C_RB_Capsule>(e)) {
            auto& cap = reg.get<cardillo::C_RB_Capsule>(e);
            const real_t newHalf = be.l * (real_t)0.5;
            if (std::abs(cap.halfLength - newHalf) > eps) {
                cap.halfLength = newHalf;
                shapeChanged = true;
            }
        }
        if (reg.any_of<cardillo::C_RB_Cylinder>(e)) {
            auto& cyl = reg.get<cardillo::C_RB_Cylinder>(e);
            const real_t newHalf = be.l * (real_t)0.5;
            if (std::abs(cyl.halfLength - newHalf) > eps) {
                cyl.halfLength = newHalf;
                shapeChanged = true;
            }
        }
        if (shapeChanged) {
            world.markStructureDirty();
        }
    }
}

void DerivedEntitySync::updateEntities(World& world, real_t dt) {
    auto& reg = world.ecs();

    Trajectory::update(world, dt);

    // Centralized rigid-state refresh for this step.
    auto poseView = reg.view<const cardillo::C_Position3>();
    for (auto e : poseView) {
        cardillo::RigidBody::updateState(reg, e);
    }

    auto view = reg.view<cardillo::C_BeamElement, const cardillo::C_Position3>();
    for (auto [e, be, pos] : view.each()) {
        (void)be;
        (void)pos;
        updateBeamElementEntity(world, e);
    }
}

}  // namespace physics
}  // namespace cardillo
