#pragma once

#include <memory>

namespace filament { class Engine; class View; }
namespace filagui { class ImGuiHelper; }
namespace frontend::camera { class CameraController; }

namespace frontend::ui {

class ImGuiLayer {
public:
    ImGuiLayer();
    ~ImGuiLayer();

    bool initialize(filament::Engine* engine, filament::View* uiView, int width, int height);
    void render(float elapsedSeconds, camera::CameraController* camera);

private:
    std::unique_ptr<filagui::ImGuiHelper> helper_;
    double fpsSmoothed_ = 0.0;
    int width_ = 0;
    int height_ = 0;
};

} // namespace filament_example::ui
