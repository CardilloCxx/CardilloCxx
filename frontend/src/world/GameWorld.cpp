#include "GameWorld.h"
#include <filament/RenderableManager.h>
#include <filament/LightManager.h>
#include <utils/EntityManager.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <cstring>

#include "geometry/Assets.h"

namespace frontend {

GameWorld::GameWorld(filament::Engine& engine, EntityPicker& entityPicker) : engine_(engine), entityPicker_(entityPicker) {
    mainScene_ = engine_.createScene();
    pickingScene_ = engine_.createScene();
    materials_ = std::make_unique<frontend::Materials>(engine_, *this);
    assets_ = std::make_unique<frontend::Assets>(engine_);
    entityPicker_.setGameWorld(this);
}

GameWorld::~GameWorld() {
    auto& em = utils::EntityManager::get();
    auto& rm = engine_.getRenderableManager();

    // Remove all tracked entities from the scenes before destroying their renderables.
    for (const auto& [id, entity] : entityToRenderable_) {
        if (!entity.isNull()) {
            mainScene_->remove(entity);
        }
    }
    for (const auto& [id, entity] : entityToPickingRenderable_) {
        if (!entity.isNull()) {
            pickingScene_->remove(entity);
        }
    }

    entityToRenderable_.clear();
    entityToPickingRenderable_.clear();

    // Destroy Lights
    auto destroyLight = [&](utils::Entity light) {
        if (!light.isNull()) {
            mainScene_->remove(light);
            em.destroy(light);
        }
    };
    destroyLight(directionalLight_);
    destroyLight(fillLight_);

    // Destroy Environment
    if (skybox_) engine_.destroy(skybox_);
    if (skyboxTexture_) engine_.destroy(skyboxTexture_);
    if (indirectLight_) engine_.destroy(indirectLight_);
    if (indirectLightTexture_) engine_.destroy(indirectLightTexture_);

    engine_.destroy(mainScene_);
    engine_.destroy(pickingScene_);
}

void GameWorld::initialize() {
    materials_->initialize();
    assets_->initialize();
    setupLights();
    setupEnvironment();
}

void GameWorld::update(double elapsedSeconds) {
    elapsedTime_ += static_cast<float>(elapsedSeconds);
//     if (indirectLight_) {
//         float rotationAngle = 0.5f * elapsedTime_;
//         
//         filament::math::mat3f rotationMatrix = filament::math::mat3f::rotation(
//             rotationAngle, 
//             filament::math::float3{0.0f, 1.0f, 0.0f} 
//         );
// 
//         indirectLight_->setRotation(rotationMatrix);
//     }
}

void GameWorld::setupLights() {
    auto& em = utils::EntityManager::get();
    
    directionalLight_ = em.create();
    filament::LightManager::Builder(filament::LightManager::Type::SUN)
        .color({1.0f, 0.88f, 0.63f})
        .intensity(3000.0f)
        .direction({-0.4f, -1.0f, -0.3f})
        .castShadows(false)
        .build(engine_, directionalLight_);
    mainScene_->addEntity(directionalLight_);

    fillLight_ = em.create();
    // filament::LightManager::Builder(filament::LightManager::Type::SUN)
    //     .color({0.7f, 0.8f, 1.0f})
    //     .intensity(3000.0f)
    //     .direction({0.2f, -1.0f, 0.2f})
    //     .castShadows(false)
    //     .build(engine_, fillLight_);
    // mainScene_->addEntity(fillLight_);
}

filament::Texture* loadTextureFromFile(filament::Engine& engine, const std::string& path) {
    std::ifstream file(path, std::ios::binary | std::ios::ate);
    if (!file) return nullptr;
    std::vector<uint8_t> data(static_cast<size_t>(file.tellg()));
    file.seekg(0);
    file.read(reinterpret_cast<char*>(data.data()), data.size());

    auto* bundle = new ktxreader::Ktx1Bundle(data.data(), static_cast<uint32_t>(data.size()));
    return ktxreader::Ktx1Reader::createTexture(&engine, bundle, false); // Ktx1Reader still wants a pointer here
}

void GameWorld::setupEnvironment() {
    // filament::Texture* skyboxTexture = loadTexture(engine_, FILAMENT_IBL_SKYBOX_PATH);
    // if (skyboxTexture) {
    //     skybox_ = filament::Skybox::Builder()
    //         .environment(skyboxTexture)
    //         .build(engine_);        
    //     mainScene_->setSkybox(skybox_);
    // }

    filament::Texture* iblTexture = loadTextureFromFile(engine_, FILAMENT_IBL_REFLECTION_PATH);
    if (iblTexture) {
        indirectLight_ = filament::IndirectLight::Builder()
            .reflections(iblTexture)
            .intensity(30000.0f)
            .build(engine_);       
        mainScene_->setIndirectLight(indirectLight_);
    }
}

utils::Entity GameWorld::createObject(PrimitiveType type, const filament::math::mat4f& transform, bool isGlass) {
    MeshData& mesh = assets_->ensureMesh(type);
    return createEntity(mesh, transform, getNextPickId(), isGlass);
}

utils::Entity GameWorld::createObject(const std::string& meshPath, const filament::math::mat4f& transform, bool isGlass) {
    MeshData& mesh = assets_->ensureMesh(meshPath);
    return createEntity(mesh, transform, getNextPickId(), isGlass);
}

utils::Entity GameWorld::createObject(const cardillo::MeshAsset& meshAsset, const filament::math::mat4f& transform, bool isGlass) {
    MeshData& mesh = assets_->ensureMesh(meshAsset);
    return createEntity(mesh, transform, getNextPickId(), isGlass);
}

utils::Entity GameWorld::createEntity(const MeshData& mesh, const filament::math::mat4f& transform, uint32_t pickId, bool isGlass) {
    auto& em = utils::EntityManager::get();
    auto& tm = engine_.getTransformManager();

    // 1. PBR Entity
    utils::Entity pbrEntity = em.create();
    entityToMeshKey_[pbrEntity.getId()] = mesh.uniqueKey;

    auto material = isGlass ? materials_->ensureGlassMaterialInstance(pbrEntity) : materials_->ensureMaterialInstance(pbrEntity);

    filament::RenderableManager::Builder(1)
        .boundingBox(mesh.bounds)
        .castShadows(true)
        .receiveShadows(true)
        .material(0, material)
        .geometry(0, mesh.primitiveType, mesh.vb, mesh.ib, 0, mesh.indexCount)
        .build(engine_, pbrEntity);

    materials_->setColor(pbrEntity, filament::math::float4(1.0f, 1.0f, 1.0f, 1.0f)); // Default to white; can be changed later
    material->setParameter("textureSwitches", filament::math::float4(0.0f, 0.0f, 0.0f, 0.0f));
    if (isGlass) { 
        materials_->setThickness(pbrEntity, mesh.bounds.getBoundingSphere().w * 0.5f);
        materials_->setDispersion(pbrEntity, 0.1f);
    } else {
        materials_->setDefaultProperties(pbrEntity); 
    }
    
    tm.create(pbrEntity);
    tm.setTransform(tm.getInstance(pbrEntity), transform);
    mainScene_->addEntity(pbrEntity);
    entityToRenderable_[pickId] = pbrEntity;


    // 2. Picking Entity
    utils::Entity pickingEntity = em.create();
    filament::RenderableManager::Builder(1)
        .boundingBox(mesh.bounds)
        .castShadows(false)
        .receiveShadows(false)
        .material(0, materials_->getPickingMaterial(pickId))
        .geometry(0, mesh.primitiveType, mesh.vb, mesh.ib, 0, mesh.indexCount)
        .build(engine_, pickingEntity);

    // Parent to PBR entity so transforms sync automatically
    tm.create(pickingEntity, tm.getInstance(pbrEntity), filament::math::mat4f());
    pickingScene_->addEntity(pickingEntity);
    entityToPickingRenderable_[pickId] = pickingEntity;

    return pbrEntity; 
}

void GameWorld::destroyInstance(utils::Entity entity) {
    const uint32_t id = entity.getId();
    auto& em = utils::EntityManager::get();
    auto& rm = engine_.getRenderableManager();

    // Destroy Main PBR Renderable
    if (auto it = entityToRenderable_.find(id); it != entityToRenderable_.end()) {
        utils::Entity r = it->second;
        mainScene_->remove(r);
        rm.destroy(r);
        em.destroy(r);
        entityToRenderable_.erase(it);
    }

    // Destroy Child Picking Renderable
    if (auto pickingIt = entityToPickingRenderable_.find(id); pickingIt != entityToPickingRenderable_.end()) {
        utils::Entity p = pickingIt->second;
        pickingScene_->remove(p);
        rm.destroy(p);
        em.destroy(p);
        entityToPickingRenderable_.erase(pickingIt);
    }
}

void GameWorld::updateTransform(utils::Entity entity, const filament::math::mat4f& transform) {
    auto& tm = engine_.getTransformManager();
    auto instance = tm.getInstance(entity);
    tm.setTransform(instance, transform);
}
}