// Contact patch convex hull of the Coal fork (patched, cmake/patches/coal-polygon-convex-hull.patch):
// duplicate, colinear and non-finite input points must not crash and give a counter-clockwise hull.
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include <coal/contact_patch/polygon_convex_hull.h>

#include <cardillo/misc/types.hpp>

using cardillo::real_t;
using Vec2 = coal::Vec2s;

namespace {

int failures = 0;

void check(bool ok, const std::string& message) {
    if (!ok) {
        ++failures;
        std::cerr << "FAILED: " << message << '\n';
    }
}

real_t cross(const Vec2& a, const Vec2& b, const Vec2& c) { return (b.x() - a.x()) * (c.y() - a.y()) - (b.y() - a.y()) * (c.x() - a.x()); }

// Strictly counter-clockwise and every input point inside or on the hull.
bool isHullOf(const std::vector<Vec2>& hull, const std::vector<Vec2>& cloud) {
    const std::size_t n = hull.size();
    if (n < 3) return false;
    for (std::size_t i = 0; i < n; ++i) {
        if (!(cross(hull[i], hull[(i + 1) % n], hull[(i + 2) % n]) > 0)) return false;
        for (const auto& p : cloud)
            if (cross(hull[i], hull[(i + 1) % n], p) < -1e-12) return false;
    }
    return true;
}

}  // namespace

int main() {
    std::vector<Vec2> hull;

    // Square with duplicated corners and colinear edge midpoints.
    const std::vector<Vec2> square = {Vec2(0, 0), Vec2(1, 0), Vec2(1, 1), Vec2(0, 1), Vec2(0, 0), Vec2(1, 1), Vec2(0.5, 0), Vec2(1, 0.5), Vec2(0.5, 0.5), Vec2(0, 0)};
    coal::computePolygonConvexHull(square, hull);
    check(hull.size() == 4 && isHullOf(hull, square), "square with duplicates and colinear points");

    // All points identical / colinear.
    const std::vector<Vec2> same(6, Vec2(0.3, -0.2));
    coal::computePolygonConvexHull(same, hull);
    check(hull.size() == 1, "identical points give one point");
    const std::vector<Vec2> line = {Vec2(0, 0), Vec2(2, 2), Vec2(1, 1), Vec2(3, 3), Vec2(1, 1)};
    coal::computePolygonConvexHull(line, hull);
    check(hull.size() == 2, "colinear points give a segment");

    // Non-finite input.
    std::vector<Vec2> nan = square;
    nan[3] = Vec2(std::numeric_limits<real_t>::quiet_NaN(), 0);
    coal::computePolygonConvexHull(nan, hull);
    check(hull.empty(), "non-finite input gives an empty hull");

    // Random clouds with many duplicates.
    std::mt19937 rng(1);
    std::uniform_int_distribution<int> grid(0, 6);
    for (int trial = 0; trial < 200; ++trial) {
        std::vector<Vec2> cloud(50);
        for (auto& p : cloud) p = Vec2(grid(rng), grid(rng));
        coal::computePolygonConvexHull(cloud, hull);
        check(isHullOf(hull, cloud), "random cloud " + std::to_string(trial));
    }

    if (failures == 0) std::cout << "test_convex_hull: all checks passed" << std::endl;
    return failures == 0 ? 0 : 1;
}
