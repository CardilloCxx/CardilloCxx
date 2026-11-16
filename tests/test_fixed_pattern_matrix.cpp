// misc::FixedPatternMatrix: the sparsity pattern holds every added entry that has been nonzero; it
// changes only when an entry becomes nonzero for the first time or entries are no longer added, and
// otherwise only the values are written.

#include <iostream>
#include <string>
#include <vector>

#include <cardillo/misc/fixed_pattern_matrix.hpp>

using namespace cardillo;

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        ++failures;
        std::cerr << "FAILED: " << message << '\n';
    }
}

struct Entry {
    int row, col;
    real_t value;
};

// One step: adds the entries to a 3 x 4 matrix; returns whether the pattern changed.
bool assemble(misc::FixedPatternMatrix& m, const std::vector<Entry>& entries) {
    m.begin(4);
    for (const auto& e : entries) m.add(e.row, e.col, e.value);
    return m.finish(3);
}

}  // namespace

int main() {
    misc::FixedPatternMatrix m;

    // First step: the nonzero entries form the pattern, the zero entry (1, 2) is left out.
    check(assemble(m, {{0, 0, 1}, {1, 2, 0}, {2, 3, 2}}), "first assembly builds the pattern");
    check(m.matrix().nonZeros() == 2 && m.matrix().coeff(0, 0) == 1 && m.matrix().coeff(2, 3) == 2, "pattern of the nonzero entries");

    // Same entries, new values: no change of the pattern.
    check(!assemble(m, {{0, 0, 3}, {1, 2, 0}, {2, 3, 4}}), "new values keep the pattern");
    check(m.matrix().coeff(0, 0) == 3 && m.matrix().coeff(2, 3) == 4, "values written");

    // (1, 2) becomes nonzero: the pattern grows; (0, 0) becomes zero and stays in it.
    check(assemble(m, {{0, 0, 0}, {1, 2, 5}, {2, 3, 4}}), "an entry that becomes nonzero grows the pattern");
    check(m.matrix().nonZeros() == 3 && m.matrix().coeff(1, 2) == 5 && m.matrix().coeff(0, 0) == 0, "zero entry stays in the pattern");
    check(!assemble(m, {{0, 0, 0}, {1, 2, 6}, {2, 3, 4}}), "zero inside the pattern changes nothing");

    // Another sequence of entries (a row added, (2, 3) no longer added): the entries still added keep
    // their place, (2, 3) leaves the pattern.
    check(assemble(m, {{0, 0, 0}, {0, 1, 7}, {1, 2, 6}}), "changed sequence of entries rebuilds the pattern");
    check(m.matrix().nonZeros() == 3 && m.matrix().coeff(0, 1) == 7 && m.matrix().coeff(2, 3) == 0, "entries no longer added leave the pattern, zeros still added stay");

    // Entries added twice at the same position are summed.
    assemble(m, {{0, 0, 1}, {0, 0, 2}, {0, 1, 7}, {1, 2, 6}});
    check(m.matrix().coeff(0, 0) == 3, "duplicate entries are summed");

    if (failures == 0) std::cout << "test_fixed_pattern_matrix: all checks passed\n";
    return failures == 0 ? 0 : 1;
}
