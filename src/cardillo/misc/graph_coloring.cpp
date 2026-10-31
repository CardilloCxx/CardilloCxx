#include "graph_coloring.hpp"

#include <algorithm>
#include <numeric>

namespace cardillo::misc {

Coloring colorGreedyWelshPowell(int numNodes, const std::vector<std::vector<int>>& adjacency) {
    Coloring result;
    result.colorOf.assign(numNodes, -1);
    if (numNodes == 0) return result;

    std::vector<int> order(numNodes);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) { return adjacency[a].size() > adjacency[b].size(); });

    // `forbiddenStamp[c] == stamp` means color c is used by a neighbor of the node currently being
    // colored. Reusing one buffer across all nodes (via a monotonically increasing stamp) avoids an
    // O(numColors) reset per node.
    std::vector<int> forbiddenStamp;
    int stamp = 0;

    for (int node : order) {
        ++stamp;
        for (int nbr : adjacency[node]) {
            int c = result.colorOf[nbr];
            if (c < 0) continue;
            if (c >= (int)forbiddenStamp.size()) forbiddenStamp.resize(c + 1, -1);
            forbiddenStamp[c] = stamp;
        }

        int chosen = 0;
        while (chosen < (int)forbiddenStamp.size() && forbiddenStamp[chosen] == stamp) ++chosen;
        result.colorOf[node] = chosen;
        result.numColors = std::max(result.numColors, chosen + 1);
    }

    result.colorClasses.resize(result.numColors);
    for (int i = 0; i < numNodes; ++i) result.colorClasses[result.colorOf[i]].push_back(i);
    return result;
}

}  // namespace cardillo::misc
