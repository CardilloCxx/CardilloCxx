#pragma once

#include <condition_variable>
#include <deque>
#include <filesystem>
#include <functional>
#include <mutex>
#include <string>
#include <thread>
#include <vector>
#include <vtkSmartPointer.h>
#include "../physics/world.hpp"
#include "mesh_generator.hpp"
#include "misc/types.hpp"
#include "pvd_writer.hpp"

class vtkPolyData;

namespace fs = std::filesystem;

namespace cardillo {
namespace collision {
struct Contact;
struct ContactManifold;
class CollisionCoal;
}
}  // namespace cardillo
namespace cardillo {
namespace physics {
class DynamicsAssembler;
}
}  // namespace cardillo

namespace cardillo::io {

// Writes simulation state as VTK PolyData (.vtp) plus a .pvd time-collection per stream, using
// the real VTK library. The in-memory vtkPolyData for each per-step frame is built synchronously
// (it has to be -- it reads live ECS/collision state), then handed off to a single background
// writer thread so the actual serialize-to-disk work doesn't block the simulation loop.
class VtkWriter {
   public:
    VtkWriter(const cardillo::config::Config& cfg);
    ~VtkWriter();

    void setOutputDir(const std::string& dir);
    void setBaseName(const std::string& name);
    void setFrequency(int freq);

    void maybeWrite(int step, real_t time, const cardillo::World& sys, cardillo::collision::CollisionCoal* collision_mgr, cardillo::misc::TimingManager* timings,
                    cardillo::physics::DynamicsAssembler* dyn = nullptr);
    void write(int step, real_t time, const cardillo::World& sys, cardillo::collision::CollisionCoal* collision_mgr, cardillo::misc::TimingManager* timings,
               cardillo::physics::DynamicsAssembler* dyn = nullptr);

    void enableContactsOutput(bool enable, const std::string& baseName) {
        m_writeContacts = enable;
        m_contactsBase = baseName;
    }
    void enableSpringsOutput(bool enable, const std::string& baseName) {
        m_writeSprings = enable;
        m_springsBase = baseName;
    }
    void enableManifoldOutput(bool enable, const std::string& baseName) {
        m_writeManifolds = enable;
        m_manifoldsBase = baseName;
    }

   private:
    const cardillo::config::Config& m_cfg;
    std::string m_outputDir;
    std::string m_baseName;
    int m_frequency{1};
    bool m_writeContacts{false};
    std::string m_contactsBase{"contacts"};
    bool m_writeSprings{false};
    std::string m_springsBase{"springs"};
    bool m_writeManifolds{false};
    std::string m_manifoldsBase{"manifolds"};
    bool m_staticGeoWritten{false};

    using EntityMesh = MeshGenerator::EntityMesh;

    // Helpers (run synchronously on the calling/sim thread)
    void enrichPressure(std::vector<EntityMesh>& meshes, const cardillo::World& sys, const std::vector<cardillo::collision::Contact>& contacts) const;
    void enrichManifolds(std::vector<cardillo::collision::ContactManifold>& manifolds, const std::vector<cardillo::collision::Contact>& contacts, real_t dt) const;
    std::string buildPath(const std::string& prefix, int step) const;
    std::string pvdPath(const std::string& prefix) const;

    vtkSmartPointer<vtkPolyData> meshesToPolyData(const std::vector<EntityMesh>& meshes) const;
    vtkSmartPointer<vtkPolyData> contactsToPolyData(const std::vector<cardillo::collision::Contact>& contacts, bool writeBodyVectors) const;
    vtkSmartPointer<vtkPolyData> springsToPolyData(const cardillo::World& sys) const;
    vtkSmartPointer<vtkPolyData> manifoldsToPolyData(const std::vector<cardillo::collision::ContactManifold>& manifolds) const;

    static void writeVtp(const vtkSmartPointer<vtkPolyData>& pd, const std::string& path);

    // Enqueues a (write .vtp + update .pvd) job for the background writer thread. `pd` is moved
    // into the job -- the caller must not keep another live reference, since VTK's own reference
    // counting is not thread-safe and this queue's mutex is what makes the ownership handoff safe.
    void enqueueFrame(vtkSmartPointer<vtkPolyData> pd, PvdWriter& pvd, const std::string& pvdFilePath, const std::string& vtpPath, const std::string& vtpBaseName, int step, real_t time);

    // Background writer thread + bounded job queue.
    std::thread m_worker;
    std::deque<std::function<void()>> m_queue;
    std::mutex m_queueMutex;
    std::condition_variable m_cvNotEmpty;
    std::condition_variable m_cvNotFull;
    bool m_stop{false};
    static constexpr std::size_t kMaxQueueDepth = 32;
    void workerLoop();

    // Per-stream .pvd state, mutated only from the (single, serial) worker thread.
    PvdWriter m_pvdGeo;
    PvdWriter m_pvdManifolds;
    PvdWriter m_pvdContacts;
    PvdWriter m_pvdSprings;
};

}  // namespace cardillo::io
