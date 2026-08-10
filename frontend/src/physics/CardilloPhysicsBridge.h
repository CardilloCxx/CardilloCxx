#pragma once

#include <entt/entt.hpp>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <physics/api/physics.hpp>
#include <math/vec3.h>
#include <filament/TransformManager.h>
#include <filament/Engine.h>

#include "scenes/SceneBase.hpp"
#include "world/GameWorld.h"
#include "world/EntityPicker.h"

namespace frontend::physics {

class CardilloPhysicsBridge {
public:
    enum class VisualKind {
        Unknown,
        Cube,
        Sphere,
        Cylinder,
        Capsule,
        Cone,
        Plane,
        Mesh,
    };

    struct RenderBinding {
        entt::entity physicsEntity{entt::null};
        utils::Entity renderEntity{};
        std::string label;
        VisualKind visualKind{VisualKind::Unknown};
        filament::math::float3 scale{1.0f, 1.0f, 1.0f};
        std::string meshPath;
    };

    explicit CardilloPhysicsBridge(GameWorld& world);
    ~CardilloPhysicsBridge() = default;

    void initialize(const std::string& configPath);
    void step(double elapsedSeconds, double maxComputationTime = 0.004);

    static cardillo::Vector3r toVector3r(const filament::math::float3& value);
    static filament::math::float3 toFloat3(const cardillo::Vector3r& value);
    static filament::math::quat toQuat(const cardillo::Quaternion4r& quaternion);
    static cardillo::Quaternion4r toQuaternion4r(const filament::math::quat& quaternion);
    static filament::math::float3 toScale3(const cardillo::Vector3r& scale);
    static filament::math::mat4f makeTransformFromPhysicsState(const cardillo::RigidBody::RigidState& state, const cardillo::Vector3r& scale);

    const cardillo::physics::PhysicsEngine& engine() const { return *engine_; }
    cardillo::physics::PhysicsEngine& engine() { return *engine_; }

private:
    void createBindingsForExistingBodies();
    void createBindingForEntity(entt::entity physicsEntity);
    void updateBindingsFromPhysics();
    void UpdateGrabInteraction();
    void updateBinding(const RenderBinding& binding);
    void updateVelocityForRenderables();

    static filament::math::mat4f makeTransformFromPhysics(const cardillo::RigidBody::RigidState& state);
    static filament::math::mat4f makeTransformFromPhysics(const cardillo::RigidBody::RigidState& state, const filament::math::float3& scale);

    void initializeScene();
    void stepSimulation(double dt);
    std::unique_ptr<SceneBase> createSceneFromConfig(const cardillo::config::Config& config);

    std::unique_ptr<SceneBase> scene_;

    GameWorld& world_;
    std::unique_ptr<cardillo::physics::PhysicsEngine> engine_;
    std::unordered_map<entt::entity, RenderBinding> bindings_;
    std::unordered_map<uint32_t, entt::entity> renderToPhysicsMap_;

    std::string configPath_;
    double accumulator_{0.0};
    double simulationTime_{0.0};
    bool initialized_{false};
};

}  // namespace filament_example::physics
