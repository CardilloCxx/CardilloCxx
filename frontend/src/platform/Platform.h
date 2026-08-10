#pragma once

#include <string>

struct GLFWwindow;

namespace frontend::platform {

struct WindowConfig {
    std::string title = "Filament Example";
    int width = 2560;
    int height = 1440;
    bool resizeable = true;
};

struct InputState {
    bool quit = false;
    bool forward = false;
    bool backward = false;
    bool left = false;
    bool right = false;
    bool sprint = false;
    bool toggleWireframe = false;
    // Mouse / pointer
    bool mouseLeftDown = false;
    bool mouseMiddleDown = false;
    bool mouseRightDown = false;
    bool shiftDown = false;
    
    double mouseDeltaX = 0.0;
    double mouseDeltaY = 0.0;
    double mouseX = 0.0;
    double mouseY = 0.0;
    double scrollDelta = 0.0;
    bool ascend = false;
    bool descend = false;
};

class Platform {
public:
    explicit Platform(WindowConfig config);
    ~Platform();

    Platform(const Platform&) = delete;
    Platform& operator=(const Platform&) = delete;

    bool initialize();
    bool shouldClose() const;
    void pollEvents();
    void swapBuffers();

    int width() const;
    int height() const;
    const InputState& inputState() const;
    InputState& inputState();
    void* nativeWindowHandle() const;

private:
    struct Impl;
    Impl* impl_;
};

} // namespace filament_example::platform
