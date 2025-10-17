#pragma once

#include <vector>
#include "types.hpp"

namespace cardillo::collision::boxbox {

// Compute OBB-OBB contacts between A and B and append to out.
// Produces 1..4 contacts for face-face, or 1 contact for edge-edge.
void collideOBB(const ObbCollider& A,
                const ObbCollider& B,
                std::vector<Contact>& out);

}
