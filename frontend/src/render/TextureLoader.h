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


namespace frontend {

class GameWorld;

class TextureLoader {

public:
    TextureLoader(filament::Engine& engine, GameWorld& gameWorld);
    ~TextureLoader();

    filament::Texture* ensureTexture(std::string& texturePath);
 
private:
    filament::Engine& engine_;
    GameWorld& gameWorld_;

    std::unordered_map<std::string, filament::Texture*> loadedTextures_;
};
}