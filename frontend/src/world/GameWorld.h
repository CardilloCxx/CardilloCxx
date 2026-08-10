#pragma once

#include <filament/Engine.h>
#include <filament/Scene.h>
#include <filament/Material.h>
#include <filament/TransformManager.h>
#include <filament/MaterialInstance.h>
#include <filament/Skybox.h>
#include <filament/IndirectLight.h>
#include <filament/Texture.h>
#include <utils/Entity.h>
#include <unordered_map>
#include <string>
#include <memory>
#include <vector>
#include <ktxreader/Ktx1Reader.h>
#include <utils/Path.h>
#include <iostream>

#include <physics/assets/assets.hpp>
#include "geometry/Geometry.h"
#include "math/ColorUtils.h"
#include "world/EntityPicker.h"
#include "render/Materials.h"
#include "geometry/Assets.h"

// Forward declare your mesh structs here
struct MeshData; 

namespace frontend {

class GameWorld {
public:
    struct frontend::EntityPicker::GrabState;

    explicit GameWorld(filament::Engine& engine, EntityPicker& entityPicker);
    ~GameWorld();

    void initialize();
    void update(double elapsedSeconds);

    filament::Scene* getMainScene() const { return mainScene_; }
    filament::Scene* getPickingScene() const { return pickingScene_; }

    utils::Entity createObject(PrimitiveType type, const filament::math::mat4f& transform, bool isGlass = false);
    utils::Entity createObject(const std::string& meshPath, const filament::math::mat4f& transform, bool isGlass = false);
    utils::Entity createObject(const cardillo::MeshAsset& mesh, const filament::math::mat4f& transform, bool isGlass = false);

    void destroyInstance(utils::Entity entity);
    void updateTransform(utils::Entity entity, const filament::math::mat4f& transform);

    utils::Entity getRenderableEntityForPickId(uint32_t pickId) const {
        auto it = entityToRenderable_.find(pickId);
        if (it != entityToRenderable_.end()) return it->second;
        return utils::Entity();
    }

    uint32_t getMeshKeyForEntity(utils::Entity entity) const {
        auto it = entityToMeshKey_.find(entity.getId());
        if (it != entityToMeshKey_.end()) return it->second;
        return 0;
    }

    bool areMaterialsReady() const { return materials_->areMaterialsReady(); }

    filament::Engine& getEngine() const { return engine_; }
    Materials& getMaterials() const { return *materials_; } 
    const frontend::EntityPicker::GrabState* getGrabState() const { return entityPicker_.getGrabState(); }

private:
    utils::Entity createEntity(const MeshData& mesh, const filament::math::mat4f& transform, uint32_t pickId, bool isGlass = false);
    uint32_t getNextPickId() { return nextPickId_++; }

    void setupLights();
    void setupEnvironment();

    filament::Engine& engine_;
    filament::Scene* mainScene_ = nullptr;
    filament::Scene* pickingScene_ = nullptr;

    // Environment & Lighting
    utils::Entity directionalLight_;
    utils::Entity fillLight_;
    filament::Skybox* skybox_ = nullptr;
    filament::IndirectLight* indirectLight_ = nullptr;
    filament::Texture* skyboxTexture_ = nullptr;
    filament::Texture* indirectLightTexture_ = nullptr;

    float elapsedTime_ = 0.0f;

    // Materials
    std::unique_ptr<frontend::Materials> materials_;
    bool scenePopulationPending_ = false;

    // Assets
    std::unique_ptr<frontend::Assets> assets_;

    // Tracking
    uint32_t nextPickId_ = 1;
    EntityPicker& entityPicker_;
    std::vector<utils::Entity> createdEntities_;
    std::unordered_map<uint32_t, uint32_t> entityToMeshKey_;
    std::unordered_map<uint32_t, utils::Entity> entityToRenderable_;
    std::unordered_map<uint32_t, utils::Entity> entityToPickingRenderable_;
};

} // namespace filament_example::scene