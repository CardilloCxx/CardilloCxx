#pragma once

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include "../../misc/types.hpp"

namespace cardillo::solver {

// Settings of nesterovIterate() (pj.max_iterations, pj.nesterov_beta_threshold, pj.nesterov_restart_limit).
struct NesterovSettings {
    int maxIterations{0};
    real_t betaThreshold{1};
    int restartLimit{0};
    int verbosity{0};      // solver.verbosity; >= 2 prints the error every 1000 sweeps
    const char* name{""};  // solver name in the log and in the divergence error
};

// Nesterov (FISTA) acceleration with adaptive restart of a fixed-point sweep, shared by projected
// Jacobi, projected Gauss-Seidel and the condensed solver. The state is the pair (a, b), which
// `sweep()` updates in place; the first `nVelocity` entries of `a` are velocities. Each iteration
// sweeps from the extrapolated point y,
//
//   x_{k+1} = sweep(y_k),   y_{k+1} = x_{k+1} + beta_{k+1} (x_{k+1} - x_k),
//
// and stops when error(v) <= 1, with v the velocities of the sweep input y_k. The momentum is reset
// (y_{k+1} = x_{k+1}) when the error grows by more than 5 %, when beta exceeds betaThreshold or
// leaves [0, 1], or when the velocity update reverses its direction, (v_y - v_{k+1}) . (v_{k+1} -
// v_k) > 0; after restartLimit resets it stays off. The buffers are rotated by swapping: besides
// the sweep, an iteration costs one pass for the extrapolation and one for the restart test. If an
// extrapolated point gives non-finite values, the iteration continues from x_k without momentum; a
// non-extrapolated one throws. Returns the number of sweeps; (a, b) hold the last sweep output.
template <class Sweep, class Error>
int nesterovIterate(VectorXr& a, VectorXr& b, Eigen::Index nVelocity, const NesterovSettings& settings, Sweep&& sweep, Error&& error) {
    VectorXr ak = a, bk = b;              // x_k
    VectorXr ay(a.size()), by(b.size());  // y_{k+1}
    VectorXr vIn(nVelocity);              // velocities of the sweep input
    double t = 1.0;
    real_t errPrev = std::numeric_limits<real_t>::infinity();
    int restarts = 0;
    bool momentum = true;
    bool extrapolated = false;  // the sweep input differs from x_k

    for (int iter = 0; iter < settings.maxIterations; ++iter) {
        vIn = a.head(nVelocity);
        sweep();
        const real_t err = error(static_cast<const VectorXr&>(vIn));
        if (settings.verbosity >= 2 && iter % 1000 == 0) std::cout << settings.name << " iter " << iter << ", error " << err << std::endl;
        if (err <= (real_t)1) return iter + 1;

        if (!std::isfinite(err) && !(a.allFinite() && b.allFinite())) {
            if (!extrapolated) throw std::runtime_error(std::string(settings.name) + " diverged");
            a = ak;
            b = bk;
            momentum = false;
            extrapolated = false;
            continue;
        }

        const double t1 = 0.5 * (1.0 + std::sqrt(4.0 * t * t + 1.0));
        const double beta = (t - 1.0) / t1;
        const auto v1 = a.head(nVelocity);
        const bool restart =
            !std::isfinite(err) || err > (real_t)1.05 * errPrev || beta > (double)settings.betaThreshold || !(beta >= 0.0 && beta <= 1.0) || (vIn - v1).dot(v1 - ak.head(nVelocity)) > 0.0;
        if (restart && ++restarts >= settings.restartLimit) momentum = false;

        if (restart || !momentum) {
            ak = a;
            bk = b;
            t = 1.0;
            extrapolated = false;
        } else {
            ay = a + (real_t)beta * (a - ak);
            by = b + (real_t)beta * (b - bk);
            ak.swap(a);  // x_k <- x_{k+1}
            bk.swap(b);
            a.swap(ay);  // sweep input <- y_{k+1}
            b.swap(by);
            t = t1;
            extrapolated = true;
        }
        errPrev = err;
    }
    if (extrapolated) {
        a.swap(ak);
        b.swap(bk);
    }
    return settings.maxIterations;
}

}  // namespace cardillo::solver
