#include "cardillo.hpp"
#include "io/vtk_writer.hpp"
#include "solver/moreau.hpp"
#include "physics/dynamics_assembler.hpp"
#include <mpi.h>
#include <Eigen/Geometry> // For AngleAxisr
// collision/collision.hpp not needed directly when using VtkWriter::maybeWriteAll

using namespace cardillo;

int main() {
    MPI_Init(nullptr, nullptr);
    int worldRank = 0, worldSize = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    PhysicsSystem sys;
    sys.setGravity(Vector3r(0,0,-9.81));

    // Basin: ground plane (z=0)
    {
        PhysicsSystem::Plane ground;
        ground.center = Vector3r(0,0,0);
        ground.normal = Vector3r(0,0,1);
        ground.up = Vector3r(0,1,0);
        ground.sizeX = 6.0; ground.sizeY = 6.0;
        sys.addObstacleBody(ground);
    }
    // const real_t basin_half = 1.0;      // basin half-extent
    // const real_t wall_height = 1.0;     // plane vertical center
    // const real_t ang_deg = 60.0;        // tilt angle from vertical
    // const real_t ang_rad = ang_deg * (M_PI / 180.0);
    // const real_t s = std::sin(ang_rad); // horizontal component magnitude
    // const real_t c = std::cos(ang_rad); // vertical component
    // const real_t hmid = 0.5 * wall_height;
    // const real_t plane_size = 3.0;
    // {
    //     PhysicsSystem::Plane px; // +X wall
    //     px.center = Vector3r( basin_half, 0, hmid);
    //     px.normal = Vector3r(-s, 0, c);
    //     px.up = Vector3r(0,1,0);
    //     px.sizeX = plane_size; px.sizeY = plane_size;
    // }
    // {
    //     PhysicsSystem::Plane nx; // -X wall
    //     nx.center = Vector3r(-basin_half, 0, hmid);
    //     nx.normal = Vector3r( s, 0, c);
    //     nx.up = Vector3r(0,1,0);
    //     nx.sizeX = plane_size; nx.sizeY = plane_size;
    // }
    // {
    //     PhysicsSystem::Plane py; // +Y wall
    //     py.center = Vector3r(0, basin_half, hmid);
    //     py.normal = Vector3r(0, -s, c);
    //     py.up = Vector3r(1,0,0);
    //     py.sizeX = plane_size; py.sizeY = plane_size;
    // }
    // {
    //     PhysicsSystem::Plane ny; // -Y wall
    //     ny.center = Vector3r(0, -basin_half, hmid);
    //     ny.normal = Vector3r(0,  s, c);
    //     ny.up = Vector3r(1,0,0);
    //     ny.sizeX = plane_size; ny.sizeY = plane_size;
    //     sys.addRigidBody(ny);
    // }

    // Billard balls in a triangular arrangement (pyramid 3D)
    const real_t r = 0.075;
    const real_t dia = 2.0 * r;
    const int rows = 20;  // works well with (1,2,3) but 4 is unstable

    const Vector3r base_center(0.0, 0.0, 0.0 + r);

    for (int k = 0; k < rows; ++k) {
        int n = rows - k;  // spheres per edge in this layer
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n - i; ++j) {
                real_t x = (2.0 * i + (( j + k))) * r;
                real_t y = std::sqrt(3.0) * (j + 1.0 / 3.0 * (k)) * r;
                real_t z = std::sqrt(6.0) * 2.0 / 3.0 * k * r;
                Vector3r pos = base_center + Vector3r(x, y, z);
                Vector3r vel = Vector3r::Zero();
                sys.addPointMass(1.0, pos, vel, r * 1.01);
            }
        }
    }
    
    // Add three walls around the base of the pyramid
    const real_t wall_thickness = 0.1;
    {
        PhysicsSystem::Cube wall1;
        wall1.center = Vector3r((rows * dia) / 2.0 - r, -wall_thickness / 2.0 - r, (dia * std::sqrt(6.0)) / 6.0);
        wall1.halfExtents = Vector3r((rows * dia) / 2.0 + r, wall_thickness / 2.0, (dia * std::sqrt(6.0)) / 3.0);
        wall1.q = Eigen::AngleAxis<real_t>(0.0, Vector3r::UnitZ());
        sys.addObstacleBody(wall1);

        PhysicsSystem::Cube wall2;
        wall2.center = Vector3r((rows * dia) / 4.0 - r, (rows * dia * std::sqrt(3.0)) / 4.0 + wall_thickness / 2.0 + r, (dia * std::sqrt(6.0)) / 6.0);
        wall2.halfExtents = Vector3r((rows * dia) / 2.0 + r, wall_thickness / 2.0, (dia * std::sqrt(6.0)) / 3.0);
        wall2.center -= Vector3r(std::cos(M_PI / 3.0) * 1.25 * r, std::sin(M_PI / 3.0) * 1.25 * r, 0.0);
        wall2.q = Eigen::AngleAxis<real_t>(M_PI / 3.0, Vector3r::UnitZ());
        sys.addObstacleBody(wall2);

        PhysicsSystem::Cube wall3;
        wall3.halfExtents = Vector3r((rows * dia) / 2.0 + r, wall_thickness / 2.0, (dia * std::sqrt(6.0)) / 3.0);
        wall3.center = Vector3r((rows * dia) * 3.0 / 4.0 - r / 2.0, (rows * dia * std::sqrt(3.0)) / 4.0 - wall_thickness / 2.0 - r, (dia * std::sqrt(6.0)) / 6.0);
        wall3.center += Vector3r(std::cos(M_PI / 3.0) * r * 1.25, std::sin(M_PI / 3.0) * r * 1.25, 0.0);
        wall3.q = Eigen::AngleAxis<real_t>(-M_PI / 3.0, Vector3r::UnitZ());
        sys.addObstacleBody(wall3);
    }
    // Bullet
    // // sys.addPointMass(1.0, Vector3r(-5.0, 0.1, ball_radius), Vector3r(50.0, -0.01, 0.0), ball_radius);
    // midpoint of pyramid: 
    sys.addPointMass(10000.0, base_center + 
        Vector3r(r * rows, rows * r / std::sqrt(3.0), 5.0),
        Vector3r(0.0, 0.0, -10.0), 10*r);

//     // Ramp example (disabled)
//     {
//         PhysicsSystem::Plane ramp;
//         ramp.center = Vector3r(3.0, 0.0, 0.5);
//         ramp.normal = Eigen::AngleAxis<real_t>(-M_PI/6, Vector3r::UnitY()) * Vector3r(0,0,1);
//         ramp.up = Eigen::AngleAxis<real_t>(-M_PI/6, Vector3r::UnitY()) * Vector3r(0,1,0);
//         ramp.sizeX = 3.0; ramp.sizeY = 3.0;
//         sys.addObstacleBody(ramp);
//     }

    // Four boundary walls around base plate
    {
        PhysicsSystem::Cube wall;
        wall.center = Vector3r(0.0, -6.0, 0.5);
        wall.halfExtents = Vector3r(6.0, 0.1, 0.5);
        wall.q = Eigen::AngleAxis<real_t>(0.0, Vector3r::UnitZ());
        sys.addObstacleBody(wall);

        wall.center = Vector3r(0.0, 6.0, 0.5);
        wall.halfExtents = Vector3r(6.0, 0.1, 0.5);
        wall.q = Eigen::AngleAxis<real_t>(0.0, Vector3r::UnitZ());
        sys.addObstacleBody(wall);

        wall.center = Vector3r(-6.0, 0.0, 0.5);
        wall.halfExtents = Vector3r(0.1, 6.0, 0.5);
        wall.q = Eigen::AngleAxis<real_t>(0.0, Vector3r::UnitZ());
        sys.addObstacleBody(wall);

        wall.center = Vector3r(6.0, 0.0, 0.5);
        wall.halfExtents = Vector3r(0.1, 6.0, 0.5);
        wall.q = Eigen::AngleAxis<real_t>(0.0, Vector3r::UnitZ());
        sys.addObstacleBody(wall);
    }

    // Writers
    // Only rank 0 writes VTK to avoid file races
    std::unique_ptr<cardillo::io::VtkWriter> writer;
    if (worldRank == 0) {
        writer = std::make_unique<cardillo::io::VtkWriter>("vtk_out", "scene", 5);
        writer->enableContactsOutput(true, "contacts");
    }

    // Simulate
    // Simulate
    const real_t dt = 1e-3;
    const real_t t1 = 5.0;
    int steps = int(t1 / dt);
    cardillo::solver::MoreauSolver solver(sys);
    real_t time = 0.0;
    if (writer) writer->maybeWrite(0, time, sys);
    for (int k = 0; k < steps; ++k) {
        solver.stepMidpoint(dt);
        time += dt;
        if (writer) writer->maybeWrite(k + 1, time, sys);
    }

    MPI_Finalize();
    return 0;
}
