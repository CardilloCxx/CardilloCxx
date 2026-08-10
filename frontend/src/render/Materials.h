#pragma once

#include <filament/Engine.h>
#include <filament/Renderer.h>
#include <filament/View.h>
#include <filament/Camera.h>
#include <filament/SwapChain.h>
#include <filament/ColorGrading.h>
#include <filament/Material.h>
#include <filament/MaterialInstance.h>

#include <vector>
#include <unordered_map>
#include <memory>

#include "TextureLoader.h"


namespace frontend {

class GameWorld;

class Materials {

public:
    enum class MaterialType {
        CONCRETE,
        METAL,
        COPPER,
        STAINLESS_STEEL,
        PLASTIC,
        WOOD,
        CARBON_FIBER,
        STONE,
        WOOD_LIGHT,
        GLASS
    };

    Materials(filament::Engine& engine, GameWorld& gameWorld);
    ~Materials();

    void initialize();
    bool areMaterialsReady() const { return materialsReady_; }

    filament::MaterialInstance* ensureMaterialInstance(utils::Entity entity);
    filament::MaterialInstance* ensureGlassMaterialInstance(utils::Entity entity);

    void setVelocity(utils::Entity entity, const filament::math::float3& velocity, const filament::math::float3& angularVelocity);
    void setColor(utils::Entity entity, const filament::math::float4& color);

    void setMetallic(utils::Entity entity, float metallic);
    void setRoughness(utils::Entity entity, float roughness);
    void setRoughnessMapScale(utils::Entity entity, float roughnessMapScale);
    void setClearCoat(utils::Entity entity, float clearCoat);
    void setClearCoatRoughness(utils::Entity entity, float clearCoatRoughness);
    void setReflectance(utils::Entity entity, float reflectance);
    void setTransmission(utils::Entity entity, float transmission);
    void setThickness(utils::Entity entity, float thickness);
    void setDispersion(utils::Entity entity, float dispersion);
    void setSpecular(utils::Entity entity, float specular);
    void setSheen(utils::Entity entity, float sheen);
    void setSheenRoughness(utils::Entity entity, float sheenRoughness);
    void setAmbientOcclusion(utils::Entity entity, float ambientOcclusion);
    void setAnisotropy(utils::Entity entity, float anisotropy);
    void setDefaultProperties(utils::Entity entity);
    
    void setUVTransform(utils::Entity entity, float offsetX, float offsetY, float scaleX, float scaleY);
    void setPBRMaps(utils::Entity entity, std::string albedoPath, std::string normalPath, std::string roughnessPath, std::string metallicPath = "", std::string aoPath = "");

    void setMaterialType(utils::Entity entity, MaterialType type, const filament::math::float4& color = filament::math::float4(1.0f, 1.0f, 1.0f, 1.0f));

    filament::MaterialInstance* getPickingMaterial(uint32_t index);

private:
    filament::Engine& engine_;
    GameWorld& gameWorld_;
    TextureLoader textureLoader_;

    filament::Texture* loadTexture(filament::Engine& engine, const std::string& path);
    filament::Material* pickingMaterial_ = nullptr;
    filament::Material* glassMaterial_ = nullptr;
    filament::Material* pbrMaterial_ = nullptr;
    std::unordered_map<uint32_t, filament::MaterialInstance*> instanceMaterialsForEntity_;

    bool materialsReady_ = false;
};
}