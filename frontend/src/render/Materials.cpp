#include "Materials.h"

#include <fstream>
#include <iostream>

#include <image/LinearImage.h>
#include <imageio/ImageDecoder.h>
#include <filament/TextureSampler.h>
#include <filament/Texture.h>
#include <algorithm>
#include <fstream>

#include "world/GameWorld.h"
#include "generated_textures.hpp"

namespace {
    std::vector<uint8_t> readMaterialPackage(const std::string& path) {
        std::ifstream input(path, std::ios::binary | std::ios::ate);
        if (!input) {
            throw std::runtime_error("Unable to open material package at: " + path);
        }
        const auto size = input.tellg();
        input.seekg(0, std::ios::beg);
        std::vector<uint8_t> bytes(static_cast<size_t>(size));
        if (size > 0) {
            input.read(reinterpret_cast<char*>(bytes.data()), size);
        }
        return bytes;
    }
}

namespace frontend {
    
Materials::Materials(filament::Engine& engine, GameWorld& gameWorld) : engine_(engine), gameWorld_(gameWorld), textureLoader_(engine, gameWorld) {}



void Materials::initialize()
{
    try {
        auto pbrData = readMaterialPackage(std::string(FILAMENT_MATERIAL_PATH) + "/pbr.filamat"); 
        pbrMaterial_ = filament::Material::Builder()
            .package(pbrData.data(), pbrData.size())
            .build(engine_);
            
        pbrMaterial_->compile(filament::Material::CompilerPriorityQueue::HIGH);
        
    } catch (const std::exception& e) {
        std::cerr << "Failed to load PBR material: " << e.what() << '\n';
    }

    try {
        auto glassData = readMaterialPackage(std::string(FILAMENT_MATERIAL_PATH) + "/glass.filamat"); 
        glassMaterial_ = filament::Material::Builder()
            .package(glassData.data(), glassData.size())
            .build(engine_);
            
        glassMaterial_->compile(filament::Material::CompilerPriorityQueue::HIGH);
        
    } catch (const std::exception& e) {
        std::cerr << "Failed to load glass material: " << e.what() << '\n';
    }

    try {
        auto pickingData = readMaterialPackage(std::string(FILAMENT_MATERIAL_PATH) + "/picking.filamat");
        pickingMaterial_ = filament::Material::Builder()
            .package(pickingData.data(), pickingData.size())
            .build(engine_);
            
        pickingMaterial_->compile(filament::Material::CompilerPriorityQueue::HIGH);
        
    } catch (const std::exception& e) {
        std::cerr << "Failed to load picking material: " << e.what() << '\n';
    }

    try { engine_.flush(); } catch (...) {}
    
    materialsReady_ = true; 
}

Materials::~Materials() {
    for (auto& pair : instanceMaterialsForEntity_) {
        if (pair.second) {
            engine_.destroy(pair.second);
        }
    }
    if (pickingMaterial_) {
        engine_.destroy(pickingMaterial_);
    }
}

filament::MaterialInstance* Materials::ensureMaterialInstance(utils::Entity entity) {

    if (instanceMaterialsForEntity_.find(entity.getId()) == instanceMaterialsForEntity_.end()) {
        auto instance = pbrMaterial_->createInstance();
        instanceMaterialsForEntity_[entity.getId()] = instance;
    }

    return instanceMaterialsForEntity_[entity.getId()];
}

filament::MaterialInstance* Materials::ensureGlassMaterialInstance(utils::Entity entity) {

    if (instanceMaterialsForEntity_.find(entity.getId()) == instanceMaterialsForEntity_.end()) {
        auto instance = glassMaterial_->createInstance();
        instanceMaterialsForEntity_[entity.getId()] = instance;
    }

    return instanceMaterialsForEntity_[entity.getId()];
}


void Materials::setVelocity(utils::Entity entity, const filament::math::float3& velocity, const filament::math::float3& angularVelocity) {
    auto matInstance = ensureMaterialInstance(entity);
    matInstance->setParameter("velocity", velocity);
    matInstance->setParameter("angularVelocity", angularVelocity);
}

void Materials::setColor(utils::Entity entity, const filament::math::float4& color) {
    auto matInstance = ensureMaterialInstance(entity);
    matInstance->setParameter("albedo", color);
}

void Materials::setMetallic(utils::Entity entity, float metallic) {
    auto matInstance = ensureMaterialInstance(entity);
    matInstance->setParameter("metallic", metallic);
}

void Materials::setRoughness(utils::Entity entity, float roughness) {
    auto matInstance = ensureMaterialInstance(entity);
    matInstance->setParameter("roughness", roughness);
}

void Materials::setRoughnessMapScale(utils::Entity entity, float roughnessMapScale) {
    auto matInstance = ensureMaterialInstance(entity);
    matInstance->setParameter("roughnessMapScale", roughnessMapScale);
}

void Materials::setClearCoat(utils::Entity entity, float clearCoat) {
    auto matInstance = ensureMaterialInstance(entity);
    matInstance->setParameter("clearCoat", clearCoat);
}

void Materials::setClearCoatRoughness(utils::Entity entity, float clearCoatRoughness) {
    auto matInstance = ensureMaterialInstance(entity);
    matInstance->setParameter("clearCoatRoughness", clearCoatRoughness);
}

void Materials::setReflectance(utils::Entity entity, float reflectance) {
    auto matInstance = ensureMaterialInstance(entity);
    matInstance->setParameter("reflectance", reflectance);
}

void Materials::setTransmission(utils::Entity entity, float transmission) {
    auto matInstance = ensureMaterialInstance(entity);
    matInstance->setParameter("transmission", transmission);
}

void Materials::setThickness(utils::Entity entity, float thickness) {
    auto matInstance = ensureMaterialInstance(entity);
    matInstance->setParameter("thickness", thickness);
}

void Materials::setDispersion(utils::Entity entity, float dispersion) {
    auto matInstance = ensureMaterialInstance(entity);
    matInstance->setParameter("dispersion", dispersion);
}

void Materials::setSpecular(utils::Entity entity, float specular) {
    auto matInstance = ensureMaterialInstance(entity);
    matInstance->setParameter("specular", specular);
}

void Materials::setSheen(utils::Entity entity, float sheen) {
    auto matInstance = ensureMaterialInstance(entity);
    matInstance->setParameter("sheen", sheen);
}

void Materials::setSheenRoughness(utils::Entity entity, float sheenRoughness) {
    auto matInstance = ensureMaterialInstance(entity);
    matInstance->setParameter("sheenRoughness", sheenRoughness);
}

void Materials::setAmbientOcclusion(utils::Entity entity, float ambientOcclusion) {
    auto matInstance = ensureMaterialInstance(entity);
    matInstance->setParameter("ambientOcclusion", ambientOcclusion);
}

void Materials::setAnisotropy(utils::Entity entity, float anisotropy) {
    auto matInstance = ensureMaterialInstance(entity);
    matInstance->setParameter("anisotropy", anisotropy);
}

void Materials::setDefaultProperties(utils::Entity entity) {
    setMetallic(entity, 0.0f);
    setRoughness(entity, 0.5f);
    setClearCoat(entity, 0.0f);
    setClearCoatRoughness(entity, 1.0f);
    setReflectance(entity, 0.5f);
    setSpecular(entity, 0.5f);
    setSheen(entity, 0.0f);
    setSheenRoughness(entity, 0.5f);
    setAmbientOcclusion(entity, 1.0f);
    setAnisotropy(entity, 0.0f);
}

void  Materials::setUVTransform(utils::Entity entity, float offsetX, float offsetY, float scaleX, float scaleY) {
    auto matInstance = ensureMaterialInstance(entity);
    matInstance->setParameter("uvCoordScale", filament::math::float4(offsetX, offsetY, scaleX, scaleY));
}

void Materials::setPBRMaps(utils::Entity entity, std::string albedoPathRel, std::string normalPathRel, std::string roughnessPathRel, std::string metallicPathRel, std::string aoPathRel) {

    std::string albedoPath = std::string(FILAMENT_TEXTURE_PATH) + "/" + albedoPathRel;
    std::string normalPath = std::string(FILAMENT_TEXTURE_PATH) + "/" + normalPathRel;
    std::string roughnessPath = std::string(FILAMENT_TEXTURE_PATH) + "/" + roughnessPathRel;
    std::string metallicPath = std::string(FILAMENT_TEXTURE_PATH) + "/" + metallicPathRel;
    std::string aoPath = std::string(FILAMENT_TEXTURE_PATH) + "/" + aoPathRel;

    bool albedoExists = std::ifstream(albedoPath).good() && !(albedoPathRel == "");
    bool normalExists = std::ifstream(normalPath).good() && !(normalPathRel == "");
    bool roughnessExists = std::ifstream(roughnessPath).good() && !(roughnessPathRel == "");
    bool metallicExists = std::ifstream(metallicPath).good() && !(metallicPathRel == "");
    bool aoExists = std::ifstream(aoPath).good() && !(aoPathRel == "");

    filament::Texture* albedoTexture = nullptr;
    filament::Texture* normalTexture = nullptr;
    filament::Texture* roughnessTexture = nullptr;
    filament::Texture* metallicTexture = nullptr;
    filament::Texture* aoTexture = nullptr;

    auto matInstance = ensureMaterialInstance(entity);
    filament::TextureSampler sampler(filament::TextureSampler::MinFilter::LINEAR, filament::TextureSampler::MagFilter::LINEAR, filament::TextureSampler::WrapMode::REPEAT);

    if (albedoExists) {
        albedoTexture = textureLoader_.ensureTexture(albedoPath);
        matInstance->setParameter("albedoMap", albedoTexture, sampler);
    } 
    if (normalExists) {
        normalTexture = textureLoader_.ensureTexture(normalPath);
        matInstance->setParameter("normalMap", normalTexture, sampler);
    } 
    if (roughnessExists) {
        roughnessTexture = textureLoader_.ensureTexture(roughnessPath);
        matInstance->setParameter("roughnessMap", roughnessTexture, sampler);
    }
    if (metallicExists) {
        metallicTexture = textureLoader_.ensureTexture(metallicPath);
        matInstance->setParameter("metallicMap", metallicTexture, sampler);
    }
    if (aoExists) {
        aoTexture = textureLoader_.ensureTexture(aoPath);
        matInstance->setParameter("aoMap", aoTexture, sampler);
    }

    matInstance->setParameter("textureSwitches", filament::math::float4(albedoExists ? 1.0f : 0.0f, normalExists ? 1.0f : 0.0f, roughnessExists ? 1.0f : 0.0f, metallicExists ? 1.0f : 0.0f));
    matInstance->setParameter("textureSwitches2", aoExists ? 1.0f : 0.0f);
}

void Materials::setMaterialType(utils::Entity entity, MaterialType type, const filament::math::float4& color) {
    auto matInstance = ensureMaterialInstance(entity);
    setColor(entity, color);

    switch (type) {
        case MaterialType::CONCRETE:
            setRoughness(entity, 0.9f);
            setPBRMaps(entity, std::string(Textures::ALBEDO_CONCRETE_KTX), std::string(Textures::NORMAL_CONCRETE_KTX), std::string(Textures::ROUGHNESS_CONCRETE_KTX));
            break;
        case MaterialType::METAL:
            setMetallic(entity, 1.0f);
            setUVTransform(entity, 0.0f, 0.0f, 16.0f, 16.0f);
            setPBRMaps(entity, "", std::string(Textures::NORMAL_SCRATCHES_KTX), std::string(Textures::ROUGHNESS_SMUDGES_KTX));
            break;
        case MaterialType::COPPER:
            setMetallic(entity, 1.0f);
            setAnisotropy(entity, 0.2f);
            setUVTransform(entity, 0.0f, 0.0f, 2.0f, 2.0f);
            setPBRMaps(entity, std::string(Textures::ALBEDO_COPPER_KTX), std::string(Textures::NORMAL_COPPER_KTX), std::string(Textures::ROUGHNESS_COPPER_KTX), std::string(Textures::METAL_COPPER_KTX));
            break;
        case MaterialType::PLASTIC:
            setClearCoat(entity, 1.0f);
            setClearCoatRoughness(entity, 0.1f);
            setPBRMaps(entity, "", std::string(Textures::NORMAL_PLASTIC_KTX), std::string(Textures::ROUGHNESS_PLASTIC_KTX));
            break;
        case MaterialType::WOOD:
            setPBRMaps(entity, std::string(Textures::ALBEDO_WOOD_KTX), std::string(Textures::NORMAL_WOOD_KTX), std::string(Textures::ROUGHNESS_WOOD_KTX));
            break;
        case MaterialType::CARBON_FIBER:
            setClearCoat(entity, 0.5f);
            setClearCoatRoughness(entity, 0.1f);
            setUVTransform(entity, 0.0f, 0.0f, 2.0f, 2.0f);
            setPBRMaps(entity, std::string(Textures::ALBEDO_CARBON_FIBER_KTX), std::string(Textures::NORMAL_CARBON_FIBER_KTX), std::string(Textures::ROUGHNESS_CARBON_FIBER_KTX));
            break;
        case MaterialType::STAINLESS_STEEL:
            setMetallic(entity, 1.0f);
            setUVTransform(entity, 0.0f, 0.0f, 2.0f, 2.0f);
            setPBRMaps(entity, std::string(Textures::ALBEDO_STAINLESS_KTX), std::string(Textures::NORMAL_STAINLESS_KTX), std::string(Textures::ROUGHNESS_STAINLESS_KTX), "", std::string(Textures::AO_STAINLESS_KTX));
            break;
        case MaterialType::STONE:
            setPBRMaps(entity, std::string(Textures::ALBEDO_STONE_KTX), std::string(Textures::NORMAL_STONE_KTX), std::string(Textures::ROUGHNESS_STONE_KTX), "", std::string(Textures::AO_STONE_KTX));
            break;
        case MaterialType::WOOD_LIGHT:
            setClearCoat(entity, 0.3f);
            setClearCoatRoughness(entity, 0.2f);
            setPBRMaps(entity, std::string(Textures::ALBEDO_WOOD_LIGHT_KTX), std::string(Textures::NORMAL_WOOD_LIGHT_KTX), std::string(Textures::ROUGHNESS_WOOD_LIGHT_KTX));
            break;        
        case MaterialType::GLASS:
            setReflectance(entity, 0.5f);
            setTransmission(entity, 1.0f);
            setDispersion(entity, 0.33f);
            setRoughnessMapScale(entity, 2.0f);
            setUVTransform(entity, 0.0f, 0.0f, 2.0f, 2.0f);
            setPBRMaps(entity, "", std::string(Textures::NORMAL_PLASTIC_KTX), std::string(Textures::ROUGHNESS_SMUDGES_KTX));
            break;
        default:
            break;
    }
}

filament::MaterialInstance* Materials::getPickingMaterial(uint32_t index) 
{
    auto instance = pickingMaterial_->createInstance(); 
    instance->setParameter("entityId", index);
    return instance;
}

}