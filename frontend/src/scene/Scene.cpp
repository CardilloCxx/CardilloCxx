#include "scene/Scene.h"

#include <filament/Engine.h>
#include <filament/Material.h>
#include <filament/IndexBuffer.h>
#include <filament/RenderableManager.h>
#include <filament/TransformManager.h>
#include <filament/VertexBuffer.h>
#include <math/vec3.h>
#include <iostream>

namespace frontend::scene {

void setupCornellBox(GameWorld& world) {
    using filament::math::mat4f;
    using filament::math::float3;
    
    // Base color constants
    const auto white = filament::math::float4(1.0, 1.0, 1.0, 1.0);
    const auto red   = filament::math::float4(1.0, 0.0, 0.0, 1.0);
    const auto green = filament::math::float4(0.0, 1.0, 0.0, 1.0);

    // --- Perfectly Aligned Box Structure ---
    // Outer box bounds are defined from -3.0 to +3.0 in X and Z, and -1.0 to +3.0 in Y.
    // Wall thickness is 0.1f. Scales represent half-extents or absolute sizes depending on your GameWorld backend.
    // Assuming standard scaling where a scale of 6.0 spans 6 units total (-3 to 3):

    // // Floor (Spans X: -3 to 3, Z: -3 to 3. Positioned slightly down so surface is at Y = -1.0)
    // const mat4f floorXform = mat4f::translation(float3(0.0f, -1.05f, 0.0f)) * mat4f::scaling(float3(6.1f, 0.1f, 6.1f));
    // auto entity = world.createObject(PrimitiveType::Cube, floorXform);
    // world.getMaterials().setColor(entity, white);

//     // Ceiling (Spans X: -3 to 3, Z: -3 to 3. Surface at Y = 3.0)
//     const mat4f ceilX = mat4f::translation(float3(0.0f, 3.05f, 0.0f)) * mat4f::scaling(float3(6.1f, 0.1f, 6.1f));
//     entity = world.createObject(PrimitiveType::Cube, ceilX);
//     world.getMaterials().setColor(entity, white);
// 
//     // Back wall (Flushed between floor, ceiling, left, and right walls)
//     const mat4f backX = mat4f::translation(float3(0.0f, 1.0f, -3.05f)) * mat4f::scaling(float3(6.1f, 4.0f, 0.1f));
//     entity = world.createObject(PrimitiveType::Cube, backX);
//     world.getMaterials().setColor(entity, white);
// 
//     // Left wall (Red - flushed edge-to-edge)
//     const mat4f leftX = mat4f::translation(float3(-3.05f, 1.0f, 0.0f)) * mat4f::scaling(float3(0.1f, 4.0f, 6.0f));
//     entity = world.createObject(PrimitiveType::Cube, leftX);
//     world.getMaterials().setColor(entity, red);
// 
//     // Right wall (Green - flushed edge-to-edge)
//     const mat4f rightX = mat4f::translation(float3(3.05f, 1.0f, 0.0f)) * mat4f::scaling(float3(0.1f, 4.0f, 6.0f));
//     entity = world.createObject(PrimitiveType::Cube, rightX);
//     world.getMaterials().setColor(entity, green);


    // --- Color-Varied Sphere Grid ---
    // Creates a 4x4 grid of spheres resting on the floor (Y = -1.0 + radius)
    const int gridCount = 4;
    const float spacing = 1.2f;
    const float startX = -1.8f;
    const float startZ = -1.8f;
    const float radius = 1.0f;

    for (int row = 0; row < gridCount; ++row) {
        for (int col = 0; col < gridCount; ++col) {
            float posX = startX + col * spacing;
            float posZ = startZ + row * spacing;
            
            uint8_t r = static_cast<uint8_t>((col * 255) / (gridCount - 1));
            uint8_t g = static_cast<uint8_t>((row * 255) / (gridCount - 1));
            uint8_t b = 200; // Constant blue baseline
            uint8_t a = 255;
            
            filament::math::float4 sphereColor = filament::math::float4(r / 255.0f, g / 255.0f, b / 255.0f, a / 255.0f);

            mat4f sphereXform = mat4f::translation(float3(posX, -1.0f + radius, posZ)) * mat4f::scaling(float3(radius));

            bool isGlass = (row == 1 && col == 2); 
            auto entity = world.createObject(PrimitiveType::Sphere, sphereXform, isGlass);
            world.getMaterials().setColor(entity, sphereColor);

            if (row ==3 && col == 3) {
                world.getMaterials().setMaterialType(entity, Materials::MaterialType::WOOD);
            }
            if (row == 3 && col == 2) {
                world.getMaterials().setMaterialType(entity, Materials::MaterialType::CARBON_FIBER);
            } 
            if (row == 3 && col == 1) {
                world.getMaterials().setMaterialType(entity, Materials::MaterialType::CONCRETE);
            }
            if (row == 3 && col == 0) {
                world.getMaterials().setMaterialType(entity, Materials::MaterialType::COPPER);
            }
            if (row == 2 && col == 0) {
                world.getMaterials().setMaterialType(entity, Materials::MaterialType::METAL);
            }
            if (row == 2 && col == 1) {
                world.getMaterials().setMaterialType(entity, Materials::MaterialType::PLASTIC, filament::math::float4(0.8f, 0.1f, 0.1f, 1.0f));
            }
            if (row == 2 && col == 2) {
                world.getMaterials().setMaterialType(entity, Materials::MaterialType::STAINLESS_STEEL);
            }
             if (row == 2 && col == 3) {
                world.getMaterials().setMaterialType(entity, Materials::MaterialType::STONE);
            }
            if (row == 1 && col == 3) {
                world.getMaterials().setMaterialType(entity, Materials::MaterialType::WOOD_LIGHT);
            }
            if (row == 1 && col == 2) {
                world.getMaterials().setMaterialType(entity, Materials::MaterialType::GLASS);
            }
        }
    }
}

} // namespace frontend::scene