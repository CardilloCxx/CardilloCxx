#include "camera/CameraController.h"

#include <filament/Camera.h>
#include <cmath>

namespace frontend::camera {

CameraController::CameraController(const frontend::platform::InputState* inputState)
    : input_(inputState) {}

void CameraController::update(double elapsedSeconds) {
    if (!input_) return;
    auto& input = const_cast<frontend::platform::InputState&>(*input_);

    // Safe helper to normalize 3D vectors
    auto normalize = [](filament::math::float3 v) {
        float len = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        if (len > 1e-6f) return filament::math::float3{v[0] / len, v[1] / len, v[2] / len};
        return v;
    };

    // Calculate base forward direction based on current orientation angles
    filament::math::float3 forward{ 
        std::cos(pitch_) * std::cos(yaw_), 
        std::sin(pitch_), 
        std::cos(pitch_) * std::sin(yaw_) 
    };
    forward = normalize(forward);

    filament::math::float3 worldUp{0.0f, 1.0f, 0.0f};
    filament::math::float3 rightVec = normalize(cross(forward, worldUp));
    filament::math::float3 upVec = normalize(cross(rightVec, forward));

    // 1. ORBIT (Right Mouse Button Drag)
    if (input.mouseRightDown) {
        const double mx = -input.mouseDeltaX;
        const double my = input.mouseDeltaY;
        
        if (mx != 0.0 || my != 0.0) {
            const float orbitSensitivity = 0.005f;
            yaw_ += static_cast<float>(mx * orbitSensitivity);
            pitch_ += static_cast<float>(my * orbitSensitivity);
            
            // Keep pitch away from gimbal lock/poles to prevent orientation flipping
            const float maxPitch = 1.49f; 
            if (pitch_ > maxPitch) pitch_ = maxPitch;
            if (pitch_ < -maxPitch) pitch_ = -maxPitch;
        }
    }

    // 2. PAN (Middle Mouse Button Drag)
    if (input.mouseMiddleDown) {
        const double mx = input.mouseDeltaX;
        const double my = input.mouseDeltaY;

        if (mx != 0.0 || my != 0.0) {
            // Scale pan speed based on distance to ensure smooth, responsive control at any scale
            const float panSensitivity = 0.0015f * distance_;
            
            // Dragging left (negative mx) should slide the camera right relative to target
            target_ -= rightVec * static_cast<float>(mx * panSensitivity);
            // Dragging down (positive my) should slide the camera up relative to target
            target_ += upVec * static_cast<float>(my * panSensitivity);
        }
    }

    // 3. ZOOM (Scroll Wheel)
    if (input.scrollDelta != 0.0) {
        const float zoomSensitivity = 0.08f;
        // Exponential zoom: zooms quickly when far, and gives fine-grain control up close
        distance_ -= static_cast<float>(input.scrollDelta) * zoomSensitivity * distance_;
        
        // Prevent camera from clipping through or reversing past the focus target
        const float minDistance = 0.1f;
        const float maxDistance = 1000.0f;
        if (distance_ < minDistance) distance_ = minDistance;
        if (distance_ > maxDistance) distance_ = maxDistance;
        
        input.scrollDelta = 0.0;
    }

    // Always clear raw mouse deltas to prevent drift/jumps once mouse buttons are released
    input.mouseDeltaX = 0.0;
    input.mouseDeltaY = 0.0;

    if (input.toggleWireframe) {
        wireframe_ = !wireframe_;
    }
}

void CameraController::applyTo(filament::Camera* camera) const {
    // Generate orientation forward vector
    const float fwdX = std::cos(pitch_) * std::cos(yaw_);
    const float fwdY = std::sin(pitch_);
    const float fwdZ = std::cos(pitch_) * std::sin(yaw_);
    const filament::math::float3 forward{fwdX, fwdY, fwdZ};

    // Camera position orbits around target_ at distance_
    const filament::math::float3 cameraPos = target_ - forward * distance_;

    // Look at target with standard world up vector
    camera->lookAt(cameraPos, target_, {0.0f, 1.0f, 0.0f});
}

} // namespace frontend::camera