#pragma once

#include <vector>

namespace cardillo::misc {

struct Coloring {
    std::vector<int> colorOf;                     // size = numNodes, colorOf[i] in [0, numColors)
    std::vector<std::vector<int>> colorClasses;    // colorClasses[c] = node indices with color c
    int numColors{0};
};

// Greedy Welsh-Powell coloring: nodes sorted by descending degree, each assigned the smallest
// color not used by any already-colored neighbor. `adjacency[i]` lists the neighbors of node i
// (must be symmetric: j in adjacency[i] iff i in adjacency[j] -- caller is responsible for this).
Coloring colorGreedyWelshPowell(int numNodes, const std::vector<std::vector<int>>& adjacency);

}  // namespace cardillo::misc
