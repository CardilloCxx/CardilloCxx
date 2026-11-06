// Compliances of the rod element (BeamConstraint) against the closed-form stiffnesses: scale factors
// of fromMaterial() enter once, torsion uses the Saint-Venant torsion constant, and the constitutive
// length of a pretensioned element is its stress-free length.
#include <cmath>
#include <iostream>
#include <string>

#include <cardillo/physics/api/physics_engine.hpp>

using namespace cardillo;
using namespace cardillo::physics;

namespace {

int failures = 0;

void check(bool ok, const std::string& message) {
    if (!ok) {
        ++failures;
        std::cerr << "FAILED: " << message << '\n';
    }
}

bool near(real_t a, real_t b) { return std::abs(a - b) <= (real_t)1e-10 * std::abs(b); }

// Row compliances (axial, shear y, shear z, torsion, bending y, bending z) of one element of length L.
VectorXr compliances(const BeamSpringParams& springs, const BeamCrossSection& section, real_t L) {
    config::Config config;
    PhysicsEngine engine(config);
    RigidProps props;
    props.mass = (real_t)1;
    const auto a = engine.addRigidBody(CubeShape(Vector3r::Constant((real_t)0.01)), RigidState(Vector3r::Zero(), Quaternion4r::Identity()), props);
    const auto b = engine.addRigidBody(CubeShape(Vector3r::Constant((real_t)0.01)), RigidState(Vector3r(L, 0, 0), Quaternion4r::Identity()), props);
    const auto index = engine.addBeamConstraint(a, b, springs, section);
    return engine.world().constraintPatterns()[index]->getConstraint().Crows;
}

}  // namespace

int main() {
    const real_t E = (real_t)2e9, nu = (real_t)0.3, G = E / (2 * (1 + nu)), L = (real_t)0.1;

    // Thin rectangle 10 x 0.5: J = b t^3 (1/3 - 0.21 r (1 - r^4/12)), r = t/b.
    const real_t b = (real_t)0.01, t = (real_t)0.0005, r = t / b;
    const auto rect = BeamCrossSection::rectangular(b, t);
    const real_t J = b * t * t * t * ((real_t)1 / 3 - (real_t)0.21 * r * (1 - r * r * r * r / 12));
    const real_t A = b * t;
    const real_t Iy = b * t * t * t / 12, Iz = b * b * b * t / 12;

    const auto scaled = BeamSpringParams::fromMaterial(E, nu, 10, 20, 3, 4, 5);
    const VectorXr c = compliances(scaled, rect, L);
    check(near(c[0], L / (10 * E * A)), "axial compliance with scale factor applied once");
    check(near(c[1], L / (20 * G * A)), "shear compliance with scale factor applied once");
    check(near(c[3], L / (3 * G * J)), "torsional compliance uses the torsion constant J");
    check(near(c[4], L / (4 * E * Iy)), "bending compliance about y");
    check(near(c[5], L / (5 * E * Iz)), "bending compliance about z");

    // Round section: J equals the polar moment.
    const real_t R = (real_t)0.002;
    const VectorXr cr = compliances(BeamSpringParams::fromMaterial(E, nu), BeamCrossSection::round(R), L);
    check(near(cr[3], L / (G * (real_t)M_PI * std::pow(R, 4) / 2)), "torsion constant of a round section");

    // Pretension: constitutive length = stress-free length L / (1 + eps).
    const real_t eps = (real_t)0.01;
    const VectorXr cp = compliances(BeamSpringParams::fromMaterial(E, nu).withAxialPretension(eps), rect, L);
    check(near(cp[0], L / (1 + eps) / (E * A)), "pretensioned element uses its stress-free length");

    // Triangle (isosceles, width w, height h): Iz = w^3 h / 48.
    const real_t w = (real_t)0.004, h = (real_t)0.003;
    check(near(BeamCrossSection::triangle(w, h).Iz(), w * w * w * h / 48), "triangle Iz about the symmetry axis");

    if (failures == 0) std::cout << "test_rod_stiffness: all checks passed" << std::endl;
    return failures == 0 ? 0 : 1;
}
