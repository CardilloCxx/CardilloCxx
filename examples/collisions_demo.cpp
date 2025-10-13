#include "cardillo.hpp"
#include "io/vtk_writer.hpp"
#include "solver/moreau.hpp"
#include <Eigen/Geometry> // For AngleAxisr
// collision/collision.hpp not needed directly when using VtkWriter::maybeWriteAll

using namespace cardillo;

int main() {
    PhysicsSystem sys;
    sys.setGravity(Vector3r(0,0,-9.81));

    // Basin: ground plane (z=0) and four inward-tilted side planes
    {
        PhysicsSystem::Plane ground;
        ground.center = Vector3r(0,0,0);
        ground.normal = Vector3r(0,0,1);
        ground.up = Vector3r(0,1,0);
        ground.sizeX = 6.0; ground.sizeY = 6.0;
        sys.addRigidBody(ground);
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
    //     sys.addRigidBody(px);
    // }
    // {
    //     PhysicsSystem::Plane nx; // -X wall
    //     nx.center = Vector3r(-basin_half, 0, hmid);
    //     nx.normal = Vector3r( s, 0, c);
    //     nx.up = Vector3r(0,1,0);
    //     nx.sizeX = plane_size; nx.sizeY = plane_size;
    //     sys.addRigidBody(nx);
    // }
    // {
    //     PhysicsSystem::Plane py; // +Y wall
    //     py.center = Vector3r(0, basin_half, hmid);
    //     py.normal = Vector3r(0, -s, c);
    //     py.up = Vector3r(1,0,0);
    //     py.sizeX = plane_size; py.sizeY = plane_size;
    //     sys.addRigidBody(py);
    // }
    // {
    //     PhysicsSystem::Plane ny; // -Y wall
    //     ny.center = Vector3r(0, -basin_half, hmid);
    //     ny.normal = Vector3r(0,  s, c);
    //     ny.up = Vector3r(1,0,0);
    //     ny.sizeX = plane_size; ny.sizeY = plane_size;
    //     sys.addRigidBody(ny);
    // }

    // Billard balls in a triangular arrangement
    const real_t ball_radius = 0.2;
    const real_t ball_dia = 2.0 * ball_radius;
    const int rows = 5;
    const real_t start_z = ball_radius + 0.01;
    const Vector3r start_pos(-ball_dia * (rows-1) / 2.0, 0, start_z);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j <= i; ++j) {
            const Vector3r pos = start_pos + Vector3r(i * ball_radius * std::sqrt(3.0), (-i + 2.0*j) * ball_radius, 0);
            const Vector3r vel = Vector3r::Zero();

            sys.addPointMass(1.0, pos + Vector3r::Random() * 0.001, vel, ball_radius);
        }
    }

    // Bullet 
    sys.addPointMass(1.0, Vector3r(-5.0, 0.1, ball_radius), Vector3r(50.0, -0.01, 0.0), ball_radius);

    // Ramp
    {
        PhysicsSystem::Plane ramp;
        ramp.center = Vector3r(3.0, 0.0, 0.5);
        ramp.normal = Eigen::AngleAxis<real_t>(-M_PI/6, Vector3r::UnitY()) * Vector3r(0,0,1);
        ramp.up = Eigen::AngleAxis<real_t>(-M_PI/6, Vector3r::UnitY()) * Vector3r(0,1,0);
        ramp.sizeX = 3.0; ramp.sizeY = 3.0;
        sys.addRigidBody(ramp);
    }

    // Walls around the base plate
    {
        PhysicsSystem::Cube wall;
        wall.center = Vector3r(0.0, -6.0, 0.5);
        wall.halfExtents = Vector3r(6.0, 0.1, 0.5);
        wall.R = Eigen::AngleAxis<real_t>(0.0, Vector3r::UnitZ()).toRotationMatrix();
        sys.addRigidBody(wall);

        wall.center = Vector3r(0.0, 6.0, 0.5);
        wall.halfExtents = Vector3r(6.0, 0.1, 0.5);
        wall.R = Eigen::AngleAxis<real_t>(0.0, Vector3r::UnitZ()).toRotationMatrix();
        sys.addRigidBody(wall);

        wall.center = Vector3r(-6.0, 0.0, 0.5);
        wall.halfExtents = Vector3r(0.1, 6.0, 0.5);
        wall.R = Eigen::AngleAxis<real_t>(0.0, Vector3r::UnitZ()).toRotationMatrix();
        sys.addRigidBody(wall);

        wall.center = Vector3r(6.0, 0.0, 0.5);
        wall.halfExtents = Vector3r(0.1, 6.0, 0.5);
        wall.R = Eigen::AngleAxis<real_t>(0.0, Vector3r::UnitZ()).toRotationMatrix();
        sys.addRigidBody(wall);
    }

    // Writers
    cardillo::io::VtkWriter writer("vtk_out", "scene", 100);
    writer.enableContactsOutput(true, "contacts");

    // Simulate
    const int steps = 6000;
    const real_t dt = 0.0005;

    cardillo::solver::MoreauSolver solver(sys);
    for (int k = 0; k < steps; ++k) {
        solver.stepMidpoint(dt);
        writer.maybeWrite(k, sys);
    }

    return 0;
}
