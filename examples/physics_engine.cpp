#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <entt/entt.hpp>
#include <iostream>
#include <utility>

/*************
 * components
 *************/
struct Mass
{
    double m;
};
struct Inertia
{
    Eigen::Matrix3d I;
};

template <int size_>
struct Coordinates
{
    static constexpr int size = size_;
    Eigen::Matrix<double, size, 1> value;
};

struct Position3 : Coordinates<3>
{
};
struct Quaternion4
{
    Eigen::Quaterniond value;
};

struct LinearVelocity3 : Coordinates<3>
{
};
struct AngularVelocity3 : Coordinates<3>
{
};

struct Force3 : Coordinates<3>
{
};
struct Moment3 : Coordinates<3>
{
};

template <int size_>
struct Index
{
    int start;
    static constexpr int size = size_;
};

struct PositionIndex3 : Index<3>
{
};
struct QuaternionIndex4 : Index<4>
{
};

struct LinearVelocityIndex3 : Index<3>
{
};
struct AngularVelocityIndex3 : Index<3>
{
};

// Tags for point masses and rigid bodies
struct PointMassTag
{
};
struct RigidBodyTag
{
};

// --- Assign DOFs ---
std::pair<int, int> assign_dofs(entt::registry& reg)
{
    int q_next = 0, v_next = 0;

    // // Note: This stores the point masses first and subsequently the regid
    // body. for (auto e : reg.view<PointMassTag>()) {
    //     reg.emplace_or_replace<PositionIndex3>(e, PositionIndex3{q_next});
    //     reg.emplace_or_replace<LinearVelocityIndex3>(e,
    //     LinearVelocityIndex3{v_next}); q_next += 3; v_next += 3;
    // }
    // for (auto e : reg.view<RigidBodyTag>()) {
    //     reg.emplace_or_replace<PositionIndex3>(e, PositionIndex3{q_next});
    //     reg.emplace_or_replace<LinearVelocityIndex3>(e,
    //     LinearVelocityIndex3{v_next}); q_next += 3; v_next += 3;
    //     reg.emplace_or_replace<QuaternionIndex4>(e,
    //     QuaternionIndex4{q_next});
    //     reg.emplace_or_replace<AngularVelocityIndex3>(e,
    //     AngularVelocityIndex3{v_next}); q_next += 4; v_next += 3;
    // }

    // Note: This stores the masses first (point masses & rigid bodies) and
    // subsequently all the inertia. Hence the mass matrix will have the
    // same structure. This might be beneficial for particle systems with a
    // few rigid bodies. Moreover, it fits well to the iteration involved in
    // the mass matrix assembly below.
    for (auto [e, pos, vel] : reg.view<Position3, LinearVelocity3>().each()) {
        reg.emplace_or_replace<PositionIndex3>(e, PositionIndex3{q_next});
        reg.emplace_or_replace<LinearVelocityIndex3>(
            e, LinearVelocityIndex3{v_next});
        q_next += 3;
        v_next += 3;
    }
    for (auto [e, quat, vel] :
         reg.view<Quaternion4, AngularVelocity3>().each()) {
        reg.emplace_or_replace<QuaternionIndex4>(e, QuaternionIndex4{q_next});
        reg.emplace_or_replace<AngularVelocityIndex3>(
            e, AngularVelocityIndex3{v_next});
        q_next += 4;
        v_next += 3;
    }
    return std::make_pair(q_next, v_next);
}

// --- Assemble Mass Matrix ---
Eigen::SparseMatrix<double> assemble_mass_matrix(const entt::registry& reg,
                                                 int total_vdofs)
{
    // TODO: Allocate the correct number of elements in order to prevent
    // resizing the vector. This is possible since we know how many masses
    // we have added.
    std::vector<Eigen::Triplet<double>> triplets;

    for (auto [e, mass, xi, vi] :
         reg.view<Mass, PositionIndex3, LinearVelocityIndex3>().each()) {
        for (int i = 0; i < 3; ++i)
            triplets.emplace_back(vi.start + i, vi.start + i, mass.m);
    }
    for (auto [e, inertia, qi, ui] :
         reg.view<Inertia, QuaternionIndex4, AngularVelocityIndex3>().each()) {
        // iterate efficient through column-major storage
        for (int j = 0; j < 3; ++j)
            for (int i = 0; i < 3; ++i)
                triplets.emplace_back(ui.start + i, ui.start + j,
                                      inertia.I(i, j));
    }
    Eigen::SparseMatrix<double> M(total_vdofs, total_vdofs);
    M.setFromTriplets(triplets.begin(), triplets.end());
    return M;
}

// --- Euler Forward Integration ---
void euler_forward_step(entt::registry& reg, double h, int total_qdofs,
                        int total_vdofs)
{
    Eigen::SparseMatrix<double> M = assemble_mass_matrix(reg, total_vdofs);
    std::cout << "M:\n" << Eigen::MatrixXd(M) << std::endl;
    Eigen::VectorXd q = Eigen::VectorXd::Zero(total_qdofs);
    Eigen::VectorXd v = Eigen::VectorXd::Zero(total_vdofs);
    Eigen::VectorXd f = Eigen::VectorXd::Zero(total_vdofs);

    // Gather v and f
    for (auto [e, vi, vel, force] :
         reg.view<LinearVelocityIndex3, LinearVelocity3, Force3>().each()) {
        v.segment<vi.size>(vi.start) = vel.value;
        f.segment<vi.size>(vi.start) = force.value;
    }
    for (auto [e, vi, vel, moment] :
         reg.view<AngularVelocityIndex3, AngularVelocity3, Moment3>().each()) {
        v.segment<vi.size>(vi.start) = vel.value;
        f.segment<vi.size>(vi.start) = moment.value;
    }

    // Solve M * dv = f
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(M);
    if (solver.info() != Eigen::Success)
        throw std::runtime_error("Mass matrix factorization failed!");
    Eigen::VectorXd dv = solver.solve(f);
    if (solver.info() != Eigen::Success)
        throw std::runtime_error("Solve failed!");

    // Update velocities and configurations
    v += h * dv;  // can be done on entity level
    for (auto [e, vi, vel, xi, pos] :
         reg.view<LinearVelocityIndex3, LinearVelocity3, PositionIndex3,
                  Position3>()
             .each()) {
        vel.value = v.segment<vi.size>(vi.start);
        pos.value += h * vel.value;
    }
    for (auto [e, ui, vel, qi, pos] :
         reg.view<AngularVelocityIndex3, AngularVelocity3, QuaternionIndex4,
                  Quaternion4>()
             .each()) {
        vel.value = v.segment<ui.size>(ui.start);
        Eigen::Quaterniond omega(0, vel.value);
        pos.value.coeffs() += 0.5 * h * (pos.value * omega).coeffs();
        pos.value.normalize();
    }
}

// --- Main ---
int main()
{
    entt::registry reg;

    auto pm = reg.create();
    reg.emplace<PointMassTag>(pm);
    reg.emplace<Mass>(pm, 1.0);
    reg.emplace<Position3>(pm, Position3{Eigen::Vector3d::Zero()});
    reg.emplace<LinearVelocity3>(pm, LinearVelocity3{Eigen::Vector3d::Zero()});
    reg.emplace<Force3>(pm, Force3{Eigen::Vector3d(0, 0, -9.81)});

    auto rb = reg.create();
    reg.emplace<RigidBodyTag>(rb);
    reg.emplace<Mass>(rb, 2.0);
    reg.emplace<Inertia>(rb, Inertia{Eigen::Matrix3d::Random()});
    reg.emplace<Position3>(rb, Position3{Eigen::Vector3d::Zero()});
    reg.emplace<LinearVelocity3>(rb, LinearVelocity3{Eigen::Vector3d::Zero()});
    reg.emplace<Quaternion4>(rb, Quaternion4{Eigen::Quaterniond::Identity()});
    reg.emplace<AngularVelocity3>(rb,
                                  AngularVelocity3{Eigen::Vector3d::Random()});
    reg.emplace<Force3>(rb, Force3{Eigen::Vector3d(0, 0, -9.81)});
    reg.emplace<Moment3>(rb, Moment3{Eigen::Vector3d(0, 0, 0)});

    auto pm2 = reg.create();
    reg.emplace<PointMassTag>(pm2);
    reg.emplace<Mass>(pm2, 3.0);
    reg.emplace<Position3>(pm2, Position3{Eigen::Vector3d::Zero()});
    reg.emplace<LinearVelocity3>(pm2, LinearVelocity3{Eigen::Vector3d::Zero()});
    reg.emplace<Force3>(pm2, Force3{Eigen::Vector3d(0, 0, -9.81)});

    auto [total_qdofs, total_vdofs] = assign_dofs(reg);

    for (int i = 0; i < 10; ++i) {
        euler_forward_step(reg, 0.01, total_qdofs, total_vdofs);
        std::cout << "Step " << i << ":\n";
        for (auto [e, x, v] :
             reg.view<PointMassTag, Position3, LinearVelocity3>().each()) {
            std::cout << "  Point mass x: " << x.value.transpose() << "\n";
            std::cout << "             v: " << v.value.transpose() << "\n";
        }
        for (auto [e, x, v, q, u] :
             reg.view<RigidBodyTag, Position3, LinearVelocity3, Quaternion4,
                      AngularVelocity3>()
                 .each()) {
            std::cout << "  Rigid body x: " << x.value.transpose() << "\n";
            std::cout << "             v: " << v.value.transpose() << "\n";
            std::cout << "             q: " << q.value.coeffs().transpose()
                      << "\n";
            std::cout << "             v: " << u.value.transpose() << "\n";
        }
    }
    return 0;
}
