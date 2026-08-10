#include "platform/Platform.h"

#include <GLFW/glfw3.h>
// Expose the X11 native accessors on Linux before including the native header.
#if defined(__linux__)
#define GLFW_EXPOSE_NATIVE_X11
#endif

#if defined(_WIN32)
#define GLFW_EXPOSE_NATIVE_WIN32
#endif

#include <GLFW/glfw3native.h>

#include <cstdint>

#include <utility>

namespace frontend::platform {

struct Platform::Impl {
    WindowConfig config;
    GLFWwindow* window = nullptr;
    InputState inputState;

    explicit Impl(WindowConfig cfg) : config(std::move(cfg)) {}
};

namespace {

void keyCallback(GLFWwindow* nativeWindow, int key, int /*scancode*/, int action, int /*mods*/) {
    auto* platform = static_cast<Platform*>(glfwGetWindowUserPointer(nativeWindow));
    if (!platform) {
        return;
    }

    auto& input = const_cast<InputState&>(platform->inputState());
    switch (key) {
        case GLFW_KEY_ESCAPE:
            input.quit = (action == GLFW_PRESS);
            break;
        case GLFW_KEY_W:
            input.forward = (action != GLFW_RELEASE);
            break;
        case GLFW_KEY_S:
            input.backward = (action != GLFW_RELEASE);
            break;
        case GLFW_KEY_A:
            input.left = (action != GLFW_RELEASE);
            break;
        case GLFW_KEY_D:
            input.right = (action != GLFW_RELEASE);
            break;
        case GLFW_KEY_LEFT_CONTROL:
            input.sprint = (action != GLFW_RELEASE);
            break;
        case GLFW_KEY_F:
            input.toggleWireframe = (action == GLFW_PRESS);
            break;
        case GLFW_KEY_SPACE:
            input.ascend = (action != GLFW_RELEASE);
            break;
        case GLFW_KEY_LEFT_SHIFT:
        case GLFW_KEY_RIGHT_SHIFT:
            input.descend = (action != GLFW_RELEASE);
            break;
        default:
            break;
    }
}

void cursorPosCallback(GLFWwindow* nativeWindow, double xpos, double ypos) {
    auto* platform = static_cast<Platform*>(glfwGetWindowUserPointer(nativeWindow));
    if (!platform) return;
    auto& input = const_cast<InputState&>(platform->inputState());
    input.mouseX = xpos;
    input.mouseY = ypos;
    static double lastX = xpos;
    static double lastY = ypos;
    input.mouseDeltaX += xpos - lastX;
    input.mouseDeltaY += ypos - lastY;
    lastX = xpos;
    lastY = ypos;
}

void mouseButtonCallback(GLFWwindow* nativeWindow, int button, int action, int /*mods*/) {
    auto* platform = static_cast<Platform*>(glfwGetWindowUserPointer(nativeWindow));
    if (!platform) return;
    auto& input = const_cast<InputState&>(platform->inputState());
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        input.mouseLeftDown = (action == GLFW_PRESS || action == GLFW_REPEAT);
    } else if (button == GLFW_MOUSE_BUTTON_MIDDLE) {
        input.mouseMiddleDown = (action == GLFW_PRESS || action == GLFW_REPEAT);
    } else if (button == GLFW_MOUSE_BUTTON_RIGHT) {
        input.mouseRightDown = (action == GLFW_PRESS || action == GLFW_REPEAT);
    }
}

void scrollCallback(GLFWwindow* nativeWindow, double xoffset, double yoffset) {
    auto* platform = static_cast<Platform*>(glfwGetWindowUserPointer(nativeWindow));
    if (!platform) return;
    auto& input = const_cast<InputState&>(platform->inputState());
    input.scrollDelta += yoffset;
}

} // namespace

Platform::Platform(WindowConfig config) : impl_(new Impl(std::move(config))) {}

Platform::~Platform() {
    if (impl_->window) {
        glfwDestroyWindow(impl_->window);
    }
    glfwTerminate();
    delete impl_;
}

bool Platform::initialize() {
    if (!glfwInit()) {
        return false;
    }

    // Use no client API: Filament will use Vulkan (or other backends) and manage rendering.
    glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API);
    glfwWindowHint(GLFW_RESIZABLE, impl_->config.resizeable ? GLFW_TRUE : GLFW_FALSE);

    impl_->window = glfwCreateWindow(impl_->config.width, impl_->config.height, impl_->config.title.c_str(), nullptr, nullptr);
    if (!impl_->window) {
        glfwTerminate();
        return false;
    }

    glfwMakeContextCurrent(impl_->window);
    glfwSwapInterval(1);

    glfwSetWindowSize(impl_->window, impl_->config.width, impl_->config.height);
    glfwSetWindowUserPointer(impl_->window, this);
    glfwSetKeyCallback(impl_->window, keyCallback);
    glfwSetCursorPosCallback(impl_->window, cursorPosCallback);
    glfwSetMouseButtonCallback(impl_->window, mouseButtonCallback);
    glfwSetScrollCallback(impl_->window, scrollCallback);

    glfwMakeContextCurrent(nullptr);

    return true;
}

bool Platform::shouldClose() const {
    return impl_->inputState.quit || glfwWindowShouldClose(impl_->window);
}

void Platform::pollEvents() {
    glfwPollEvents();
}

void Platform::swapBuffers() {
    if (glfwGetCurrentContext()) {
        glfwSwapBuffers(impl_->window);
    }
}

int Platform::width() const {
    int width = 0;
    glfwGetWindowSize(impl_->window, &width, nullptr);
    return width;
}

int Platform::height() const {
    int height = 0;
    glfwGetWindowSize(impl_->window, nullptr, &height);
    return height;
}

const InputState& Platform::inputState() const {
    return impl_->inputState;
}

InputState& Platform::inputState() {
    return impl_->inputState;
}

void* Platform::nativeWindowHandle() const {
#if defined(__linux__)
    return reinterpret_cast<void*>(static_cast<uintptr_t>(glfwGetX11Window(impl_->window)));
#elif defined(_WIN32)
    return static_cast<void*>(glfwGetWin32Window(impl_->window));
#else
    return nullptr;
#endif
}
} // namespace filament_example::platform
