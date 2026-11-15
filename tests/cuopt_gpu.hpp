#pragma once

// cuOpt solves only on an NVIDIA GPU. Without one, the tests skip their cuOpt checks; with
// CARDILLO_REQUIRE_GPU=1 a missing GPU is an error, so that a broken GPU runner cannot pass.

#include <cstdlib>
#include <iostream>
#include <string>

#include <cardillo/physics/solver/conic/cuopt/cuopt_backend.hpp>

namespace cuopt_test {

inline bool gpuRequired() {
    const char* v = std::getenv("CARDILLO_REQUIRE_GPU");
    return v != nullptr && std::string(v) == "1";
}

/// True if cuOpt can solve on this machine; otherwise prints why and whether the caller skips or fails.
inline bool gpuAvailable(const std::string& test) {
    std::string reason;
    if (cardillo::solver::cuopt::deviceAvailable(&reason)) return true;
    std::cout << test << ": cuOpt cannot solve on this machine (" << reason << "); " << (gpuRequired() ? "CARDILLO_REQUIRE_GPU=1, failing" : "skipping the cuOpt checks") << std::endl;
    return false;
}

}  // namespace cuopt_test
