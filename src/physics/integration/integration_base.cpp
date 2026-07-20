#include "integration_base.hpp"

namespace cardillo {
namespace integration {

void IntegrationBase::explicitPositionUpdate(World& world, real_t h) {
    auto _sc = m_timings.scope(misc::TimingManager::TimerId::Integration);
    auto& reg = world.ecs();

    auto position_view = reg.view<C_Position3, const C_LinearVelocity3>();
    for (auto [e, pos, vel] : position_view.each()) {
        (void)e;
        pos.value += h * vel.value;
    }

    auto orientation_view = reg.view<C_Orientation, const C_AngularVelocity3>();
    for (auto [e, orientation, angularVel] : orientation_view.each()) {
        (void)e;
        const Vector3r& omega = angularVel.value;
        const Quaternion4r q_prev = orientation.value;

        Vector4r& P = orientation.value.coeffs();
        real_t w = P(3);
        P(3) -= h * 0.5 * P.head<3>().dot(omega);
        P.head<3>() += h * 0.5 * (w * omega + P.head<3>().cross(omega));

        orientation.value = MathHelper::alignQuaternionTo(orientation.value, q_prev);
    }
}

void IntegrationBase::linearImplicitPositionUpdate(World& world, real_t h) {
    auto _sc = m_timings.scope(misc::TimingManager::TimerId::Integration);
    auto& reg = world.ecs();

    auto position_view = reg.view<C_Position3, const C_LinearVelocity3>();
    for (auto [e, pos, vel] : position_view.each()) {
        (void)e;
        pos.value += h * vel.value;
    }

    auto orientation_view = reg.view<C_Orientation, const C_AngularVelocity3>();
    for (auto [e, orientation, angularVel] : orientation_view.each()) {
        (void)e;
        const Vector3r& omega = angularVel.value;
        const Quaternion4r q_prev = orientation.value;

        Matrix44r D = Matrix44r::Zero();
        D.topLeftCorner<3, 3>() = SkewSymmetricMatrix3r(-omega);
        D.topRightCorner<3, 1>() = omega;
        D.bottomLeftCorner<1, 3>() = -omega.transpose();

        const Matrix44r A = Matrix44r::Identity() - 0.5 * h * D;

        Vector4r& P = orientation.value.coeffs();
        P = A.inverse() * P;

        orientation.value = MathHelper::alignQuaternionTo(orientation.value, q_prev);
    }
}

}  // namespace integration
}  // namespace cardillo
