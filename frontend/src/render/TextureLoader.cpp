#include "TextureLoader.h"
#include <filament/Engine.h>
#include <filament/Texture.h>
#include <ktxreader/Ktx1Reader.h>
#include <fstream>
#include <iostream>

namespace frontend {
    TextureLoader::TextureLoader(filament::Engine& engine, GameWorld& gameWorld) : engine_(engine), gameWorld_(gameWorld) {}

    TextureLoader::~TextureLoader() {
        for (auto& [path, texture] : loadedTextures_) {
            if (texture) engine_.destroy(texture);
        }
    }

    filament::Texture* loadTexture(filament::Engine& engine, const std::string& path) {
        std::ifstream file(path, std::ios::binary | std::ios::ate);
        if (!file) return nullptr;
        std::vector<uint8_t> data(static_cast<size_t>(file.tellg()));
        file.seekg(0);
        file.read(reinterpret_cast<char*>(data.data()), data.size());

        auto* bundle = new ktxreader::Ktx1Bundle(data.data(), static_cast<uint32_t>(data.size()));
        return ktxreader::Ktx1Reader::createTexture(&engine, bundle, false);
    }

    filament::Texture* TextureLoader::ensureTexture(std::string& texturePath) {
        auto it = loadedTextures_.find(texturePath);
        if (it != loadedTextures_.end()) return it->second;

        filament::Texture* texture = loadTexture(engine_, texturePath);
        if (!texture) {
            std::cerr << "Failed to load texture: " << texturePath << std::endl;
            return nullptr;
        }

        loadedTextures_[texturePath] = texture;
        return texture;
    }

}