#pragma once

#include "platform/Platform.h"
#include <math/vec3.h>

namespace filament { class Camera; }

namespace frontend::camera {

class CameraController {
public:
    explicit CameraController(const platform::InputState* inputState);
    void update(double elapsedSeconds);
    void applyTo(filament::Camera* camera) const;

    float movementSpeed() const { return movementSpeed_; }
    void setMovementSpeed(float s) { movementSpeed_ = s; }
    bool wireframe() const { return wireframe_; }
    void setWireframe(bool w) { wireframe_ = w; }

private:
    const frontend::platform::InputState* input_ = nullptr;
    
    filament::math::float3 target_{0.0f, 0.0f, 0.0f};   // The point the camera orbits around
    float distance_{5.0f};                             // Current distance from the target
    float yaw_{0.0f};                                  // Horizontal rotation angle
    float pitch_{0.0f};                                 // Vertical rotation angle
    bool wireframe_ = false;

    filament::math::float3 position_{-1.0f, 0.0f, 0.0f};
    float movementSpeed_ = 0.5f;
};

} // namespace frontend::camera
