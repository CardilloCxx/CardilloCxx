#include "physics/CardilloPhysicsBridge.h"

#include <Eigen/Geometry>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

#include <physics/api/physics.hpp>
#include <config/path.hpp>

#include "scenes/wilberforce/WilberforcePendulum.hpp"
#include "scenes/cradle/NewtonsCradleScene.hpp"
#include "scenes/domino/DominoScene.hpp"
#include "scenes/woodpecker/WoodpeckerScene.hpp"
#include "scenes/cantilever/CantileverScene.hpp"
#include "scenes/billiard/BilliardScene.hpp"
#include "scenes/cable/CableScene.hpp"
#include "scenes/jenga/scene.hpp"

#include "geometry/Assets.h"

namespace frontend::physics {
namespace {

const char* kDefaultSceneConfig = "examples/scenes/cradle/scene.config"; 
//"examples/scenes/wilberforce/scene.config";

std::string resolveConfigPath(const std::string& path) {
    const std::filesystem::path candidate(path);
    if (!path.empty() && candidate.is_absolute() && std::filesystem::exists(candidate)) {
        return candidate.string();
    }

    const std::filesystem::path projectCandidate = std::filesystem::path(PROJECT_SOURCE_DIR) / candidate;
    if (!path.empty() && std::filesystem::exists(projectCandidate)) {
        return projectCandidate.string();
    }

    return path;
}

filament::math::float3 CardilloPhysicsBridge::toFloat3(const cardillo::Vector3r& value) {
    return filament::math::float3(
        static_cast<float>(value.x()),
        static_cast<float>(value.z()),     
        static_cast<float>(-value.y())     
    );
}

filament::math::quat CardilloPhysicsBridge::toQuat(const cardillo::Quaternion4r& quaternion) {
    return filament::math::quat(
        static_cast<float>(quaternion.w()), // W remains unchanged
        static_cast<float>(quaternion.x()), 
        static_cast<float>(quaternion.z()), 
        static_cast<float>(-quaternion.y()) 
    );
}

cardillo::Vector3r CardilloPhysicsBridge::toVector3r(const filament::math::float3& value) {
    return cardillo::Vector3r(
        static_cast<double>(value.x),
        static_cast<double>(-value.z), // Reversing the negation
        static_cast<double>(value.y)   // Swapping z and y back
    );
}

cardillo::Quaternion4r CardilloPhysicsBridge::toQuaternion4r(const filament::math::quat& quaternion) {
    return cardillo::Quaternion4r(
        static_cast<double>(quaternion.w),
        static_cast<double>(quaternion.x),
        static_cast<double>(-quaternion.z), // Reversing the negation
        static_cast<double>(quaternion.y)   // Swapping z and y back
    );
}

filament::math::float3 CardilloPhysicsBridge::toScale3(const cardillo::Vector3r& scale) {
    return filament::math::float3(
        static_cast<float>(scale.x()),
        static_cast<float>(scale.y()),
        static_cast<float>(scale.z())
    );
}

filament::math::mat4f CardilloPhysicsBridge::makeTransformFromPhysicsState(const cardillo::RigidBody::RigidState& state, const cardillo::Vector3r& scale) {
    const auto translation = toFloat3(state.position);
    const auto rotation = toQuat(state.orientation);
    const auto scaling = toScale3(scale);

    return filament::math::mat4f::translation(translation) *
           filament::math::mat4f(rotation) *
           filament::math::mat4f::scaling(scaling);
}

std::string extractLabel(const entt::registry& registry, entt::entity entity) {
    if (entity == entt::null || !registry.valid(entity)) {
        return {};
    }

    if (registry.any_of<cardillo::C_TrackTag>(entity)) {
        return registry.get<cardillo::C_TrackTag>(entity).name;
    }

    return std::string("entity_") + std::to_string(static_cast<uint32_t>(entity));
}

}  // namespace

CardilloPhysicsBridge::CardilloPhysicsBridge(GameWorld& world) : world_(world) {}

void CardilloPhysicsBridge::initialize(const std::string& configPath) {
    configPath_ = configPath.empty() ? kDefaultSceneConfig : configPath;
    const std::string resolvedPath = resolveConfigPath(configPath_);

    auto config = cardillo::config::ConfigReader::fromFile(resolvedPath);
    config.sim_T = -1.0;
    config.output_interval_steps = 100000000;
    engine_ = std::make_unique<cardillo::physics::PhysicsEngine>(config);

    initializeScene();
    std::cout << "Physics bridge initialized from " << configPath_ << std::endl;

    initialized_ = true;
}
void CardilloPhysicsBridge::UpdateGrabInteraction() {
    const EntityPicker::GrabState* grabState = world_.getGrabState();

    if (!grabState || !grabState->isGrabbing || grabState->hoveredEntityId == 0) {
        return;
    }

    auto it = renderToPhysicsMap_.find(grabState->hoveredEntityRenderableId);
    if (it == renderToPhysicsMap_.end()) {
        return;
    }

    entt::entity hoveredPhysicsEntity = it->second;

    if (engine_->isStatic(hoveredPhysicsEntity)) {
        return;
    }

    const double entityMass = engine_->getMass(hoveredPhysicsEntity)(0, 0);
    const double stiffness = 300.0 * entityMass;
    const double damping = 30.0 * entityMass;

    filament::TransformManager& transformManager = world_.getEngine().getTransformManager();
    utils::Entity entity = world_.getRenderableEntityForPickId(grabState->hoveredEntityId);
    auto transformInstance = transformManager.getInstance(entity);
    
    const filament::math::mat4f worldTransform = transformManager.getWorldTransform(transformInstance);
    const filament::math::float3 entityPosition = worldTransform[3].xyz;
    const filament::math::mat3f entityRotation = worldTransform.upperLeft();

    filament::math::mat3f pureRotation(
        normalize(entityRotation[0]),
        normalize(entityRotation[1]),
        normalize(entityRotation[2])
    );

    filament::math::float3 grabPointInWorld = entityPosition + (pureRotation * grabState->inBodyGrabPoint);

    cardillo::Vector3r linearVelWorld = engine_->getLinearVelocity(hoveredPhysicsEntity);
    cardillo::Vector3r angularVelBody = engine_->getAngularVelocity(hoveredPhysicsEntity);

    cardillo::Vector3r r_local_physics = toVector3r(grabState->inBodyGrabPoint);
    cardillo::Vector3r tangentialVelBody = angularVelBody.cross(r_local_physics);

    filament::math::float3 tangVelBodyFilament = toFloat3(tangentialVelBody);
    filament::math::float3 tangVelWorldFilament = pureRotation * tangVelBodyFilament;
    cardillo::Vector3r tangentialVelWorld = toVector3r(tangVelWorldFilament);

    cardillo::Vector3r pointVelWorld = linearVelWorld + tangentialVelWorld;     

    cardillo::Vector3r currentWorldGrabPoint = toVector3r(grabPointInWorld);
    cardillo::Vector3r targetWorldPos = toVector3r(grabState->mouseWorldPos);

    cardillo::Vector3r displacement = targetWorldPos - currentWorldGrabPoint;
    cardillo::Vector3r springForce = (displacement * stiffness) - (pointVelWorld * damping);

    engine_->applyForceAt(hoveredPhysicsEntity, springForce, r_local_physics);
}

void CardilloPhysicsBridge::step(double elapsedSeconds, double maxComputationTime) {
    if (!initialized_ || !engine_) {
        return;
    }

    accumulator_ += elapsedSeconds;
    const double dt = engine_->world().config().sim_dt;

    if (dt <= 0.0) {
        return;
    }

    const auto startTime = std::chrono::steady_clock::now();

    while (accumulator_ >= dt) {
        const auto currentTime = std::chrono::steady_clock::now();
        const double elapsedTime = std::chrono::duration<double>(currentTime - startTime).count();

        if (elapsedTime >= maxComputationTime) {
            break;
        }

        stepSimulation(dt);
        accumulator_ -= dt;
    }

    updateBindingsFromPhysics();
    updateVelocityForRenderables();
}

void CardilloPhysicsBridge::createBindingsForExistingBodies() {
    auto& registry = engine_->ecs();
    auto view = registry.view<cardillo::C_Position3, cardillo::C_Orientation>();
    view.each([this, &registry](entt::entity entity, const cardillo::C_Position3&, const cardillo::C_Orientation&) {
        if (registry.any_of<cardillo::C_VisualObject>(entity)) {
            createBindingForEntity(entity);
        }
    });
}

void CardilloPhysicsBridge::createBindingForEntity(entt::entity physicsEntity) {
    if (physicsEntity == entt::null || !engine_->ecs().valid(physicsEntity)) {
        return;
    }

    auto& registry = engine_->ecs();
    if (bindings_.contains(physicsEntity)) {
        return;
    }

    const auto state = cardillo::RigidBody::getState(registry, physicsEntity);
    const auto label = extractLabel(registry, physicsEntity);

    RenderBinding binding{physicsEntity, utils::Entity{}, label};
    binding.scale = filament::math::float3{1.0f, 1.0f, 1.0f};
// 
    const float colorScalar = (float) physicsEntity / (54.0f);
    
    const auto transferFunc = [](float value) -> filament::math::float3 {
        const float r = std::sin(value * 6.28318530718f + 0.0f) * 0.5f + 0.5f;
        const float g = std::sin(value * 6.28318530718f + 2.09439510239f) * 0.5f + 0.5f;
        const float b = std::sin(value * 6.28318530718f + 4.18879020479f) * 0.5f + 0.5f;

        return filament::math::float3{r, g, b};
    }; 

    auto colorRgb = transferFunc(colorScalar);

    if (cardillo::RigidBody::isStatic(engine_->ecs(), physicsEntity)) {
        colorRgb = filament::math::float3{0.5f, 0.5f, 0.5f};
    }

    if (registry.any_of<cardillo::C_MeshVisualTag>(physicsEntity) && registry.any_of<cardillo::C_Mesh>(physicsEntity)) {
        const auto& asset = engine_->world().getMeshAsset(physicsEntity);
        binding.visualKind = VisualKind::Mesh;
        binding.renderEntity = world_.createObject(asset, makeTransformFromPhysicsState(state, cardillo::Vector3r(1.0, 1.0, 1.0)));
    } else if (registry.any_of<cardillo::C_CylinderVisualTag>(physicsEntity)) {
        binding.visualKind = VisualKind::Cylinder;
        if (registry.any_of<cardillo::C_Cylinder>(physicsEntity)) {
            const auto& cylinder = registry.get<cardillo::C_Cylinder>(physicsEntity);
            binding.scale = 2.0 * filament::math::float3{static_cast<float>(cylinder.radius), static_cast<float>(cylinder.halfLength), static_cast<float>(cylinder.radius)};
        }
        binding.renderEntity = world_.createObject(PrimitiveType::Cylinder, makeTransformFromPhysicsState(state, cardillo::Vector3r(1.0, 1.0, 1.0)));
    } else if (registry.any_of<cardillo::C_CapsuleVisualTag>(physicsEntity)) {
        binding.visualKind = VisualKind::Cylinder;
        if (registry.any_of<cardillo::C_Capsule>(physicsEntity)) {
            const auto& capsule = registry.get<cardillo::C_Capsule>(physicsEntity);
            binding.scale = 2.0 * filament::math::float3{static_cast<float>(capsule.radius), static_cast<float>(capsule.halfLength), static_cast<float>(capsule.radius)};
        }
        binding.renderEntity = world_.createObject(PrimitiveType::Cylinder, makeTransformFromPhysicsState(state, cardillo::Vector3r(1.0, 1.0, 1.0)));
    } else if (registry.any_of<cardillo::C_ConeVisualTag>(physicsEntity)) {
        binding.visualKind = VisualKind::Cone;
        if (registry.any_of<cardillo::C_Cone>(physicsEntity)) {
            const auto& cone = registry.get<cardillo::C_Cone>(physicsEntity);
            binding.scale = 2.0 * filament::math::float3{static_cast<float>(cone.radius),  static_cast<float>(cone.height), static_cast<float>(cone.radius) };
        }
        binding.renderEntity = world_.createObject(PrimitiveType::Cone, makeTransformFromPhysicsState(state, cardillo::Vector3r(1.0, 1.0, 1.0)));
    } else if (registry.any_of<cardillo::C_PointVisualTag>(physicsEntity)) {
        binding.visualKind = VisualKind::Sphere;
        if (registry.any_of<cardillo::C_Radius>(physicsEntity)) {
            const auto& radius = registry.get<cardillo::C_Radius>(physicsEntity);
            binding.scale = 2.0 * filament::math::float3{static_cast<float>(radius.r), static_cast<float>(radius.r), static_cast<float>(radius.r)};
        }
        binding.renderEntity = world_.createObject(PrimitiveType::Sphere, makeTransformFromPhysicsState(state, cardillo::Vector3r(1.0, 1.0, 1.0)));
    } else if (registry.any_of<cardillo::C_CubeVisualTag>(physicsEntity)) {
        binding.visualKind = VisualKind::Cube;
        if (registry.any_of<cardillo::C_Cube>(physicsEntity)) {
            const auto& cube = registry.get<cardillo::C_Cube>(physicsEntity);
            binding.scale = filament::math::float3{static_cast<float>(cube.halfExtents.x() * 2.0), static_cast<float>(cube.halfExtents.z() * 2.0), static_cast<float>(cube.halfExtents.y() * 2.0)};
        }
        binding.renderEntity = world_.createObject(PrimitiveType::Cube, makeTransformFromPhysicsState(state, cardillo::Vector3r(1.0, 1.0, 1.0)));
    } else {
        binding.visualKind = VisualKind::Cube;
        binding.renderEntity = world_.createObject(PrimitiveType::Cube, makeTransformFromPhysicsState(state, cardillo::Vector3r(0.01, 0.01, 0.01)));
    }

    Materials::MaterialType materialType = ((int)physicsEntity == 0 ? Materials::MaterialType::WOOD : Materials::MaterialType::PLASTIC);

    world_.getMaterials().setMaterialType(binding.renderEntity, materialType);
    float uvScale = 2.0;
    world_.getMaterials().setUVTransform(binding.renderEntity, 0.0f, 0.0f, uvScale * binding.scale.x, uvScale * binding.scale.z);
    world_.getMaterials().setColor(binding.renderEntity, filament::math::float4(colorRgb.x, colorRgb.y, colorRgb.z, 1.0f));

    bindings_[physicsEntity] = std::move(binding);
    renderToPhysicsMap_[binding.renderEntity.getId()] = physicsEntity;
}

void CardilloPhysicsBridge::updateBindingsFromPhysics() {
    for (const auto& [physicsEntity, binding] : bindings_) {
        updateBinding(binding);
    }
}

void CardilloPhysicsBridge::updateBinding(const RenderBinding& binding) {
    if (binding.physicsEntity == entt::null || !engine_->ecs().valid(binding.physicsEntity)) {
        return;
    }

    const auto state = cardillo::RigidBody::getState(engine_->ecs(), binding.physicsEntity);
    const auto transform = makeTransformFromPhysicsState(state, cardillo::Vector3r(binding.scale.x, binding.scale.y, binding.scale.z));
    world_.updateTransform(binding.renderEntity, transform);
}

void CardilloPhysicsBridge::updateVelocityForRenderables() {
    for (const auto& [physicsEntity, binding] : bindings_) {
        if (binding.physicsEntity == entt::null || !engine_->ecs().valid(binding.physicsEntity)) {
            continue;
        }

        const auto state = cardillo::RigidBody::getState(engine_->ecs(), binding.physicsEntity);
        world_.getMaterials().setVelocity(binding.renderEntity, toFloat3(state.linearVelocity), toFloat3(state.angularVelocity));
    }
}

filament::math::mat4f CardilloPhysicsBridge::makeTransformFromPhysics(const cardillo::RigidBody::RigidState& state) {
    return makeTransformFromPhysicsState(state, cardillo::Vector3r(1.0, 1.0, 1.0));
}

filament::math::mat4f CardilloPhysicsBridge::makeTransformFromPhysics(const cardillo::RigidBody::RigidState& state, const filament::math::float3& scale) {
    return makeTransformFromPhysicsState(state, cardillo::Vector3r(scale.x, scale.y, scale.z));
}

void CardilloPhysicsBridge::initializeScene() {
    scene_ = createSceneFromConfig(engine_->world().config());
    if (!scene_) {
        throw std::runtime_error("No Cardillo scene could be instantiated for the current configuration.");
    }

    scene_->populate(*engine_);
    createBindingsForExistingBodies();
}

void CardilloPhysicsBridge::stepSimulation(double dt) {
    if (!scene_) {
        return;
    }

    scene_->updateScene(*engine_, simulationTime_, dt);
    UpdateGrabInteraction();
    engine_->step();
    updateBindingsFromPhysics();
    simulationTime_ += dt;
}

std::unique_ptr<SceneBase> CardilloPhysicsBridge::createSceneFromConfig(const cardillo::config::Config& config) {
    if (config.scene_name == "wilberforce") {
        return std::make_unique<WilberforcePendulumScene>();
    }else if (config.scene_name == "cradle") {
        return std::make_unique<NewtonsCradleScene>();
    }else if (config.scene_name == "domino") {
        return std::make_unique<DominoScene>();
    }else if (config.scene_name == "woodpecker") {
        return std::make_unique<WoodpeckerScene>();
    }else if (config.scene_name == "cantilever") {
        return std::make_unique<CantileverScene>();
    }else if (config.scene_name == "billiard") {
        return std::make_unique<BilliardScene>();
    }else if (config.scene_name == "cable") {
        return std::make_unique<CableScene>();
    }else if (config.scene_name == "jenga") {
        return std::make_unique<JengaScene>();
    }

    throw std::runtime_error("Unsupported Cardillo scene: " + config.scene_name);
}

}  // namespace frontend::physics
