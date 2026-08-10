#include "ui/ImGuiLayer.h"

#include <filagui/ImGuiHelper.h>
#include <imgui.h>
#include <filament/Engine.h>
#include <filament/View.h>
#include "camera/CameraController.h"
#include <utils/Path.h>

namespace frontend::ui {

ImGuiLayer::ImGuiLayer() = default;
ImGuiLayer::~ImGuiLayer() = default;

bool ImGuiLayer::initialize(filament::Engine* engine, filament::View* uiView, int width, int height) {
    if (!engine || !uiView) return false;
    helper_ = std::make_unique<filagui::ImGuiHelper>(engine, uiView, utils::Path(""));
    helper_->setDisplaySize(width, height);
    width_ = width; height_ = height;
    return true;
}

void ImGuiLayer::render(float elapsedSeconds, camera::CameraController* camera) {
    if (!helper_) return;
    helper_->render(elapsedSeconds, [this, camera, elapsedSeconds](filament::Engine* engine, filament::View* view){
        ImGui::SetNextWindowPos(ImVec2(10,10), ImGuiCond_Always);
        ImGui::Begin("Settings", nullptr, ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoSavedSettings);
        if (camera) {
            bool wf = camera->wireframe();
            if (ImGui::Checkbox("Wireframe", &wf)) camera->setWireframe(wf);
            float ms = camera->movementSpeed();
            if (ImGui::SliderFloat("Move Speed", &ms, 0.1f, 20.0f)) camera->setMovementSpeed(ms);
        }
        ImGui::End();

        // FPS overlay top-right
        double fps = elapsedSeconds > 0.0 ? (1.0 / elapsedSeconds) : 0.0;
        if (fpsSmoothed_ == 0.0) fpsSmoothed_ = fps; else fpsSmoothed_ = fpsSmoothed_ * 0.95 + fps * 0.05;

        ImGui::SetNextWindowBgAlpha(0.3f);
        ImGuiWindowFlags fpsFlags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings;
        ImGui::SetNextWindowPos(ImVec2((float)width_ - 10.0f, 10.0f), ImGuiCond_Always, ImVec2(1.0f, 0.0f));
        ImGui::Begin("##fps", nullptr, fpsFlags);
        ImGui::Text("FPS: %.1f", fpsSmoothed_);
        ImGui::End();
    });
}

} // namespace filament_example::ui
