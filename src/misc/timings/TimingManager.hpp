// TimingManager.hpp - lightweight timing aggregation for subsystems
#pragma once
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace cardillo {
namespace misc {

class TimingManager {
   public:
    enum class TimerId {
        CollisionBroadphase,
        CollisionMatching,
        CollisionNarrowphase,
        UpdateEntities,
        BuildAndFactorS,
        RebuildConstraintJacobians,
        RebuildContactJacobians,
        Integration,
        ProjectedJacobi,
        ProjectedJacobiSweep,
        ProjectedNewton,
        ProjectedNewtonSetup,
        QocoSolve,
        ClarabelSolve,
        CollisionMakeContact,
        DisableCollisionPairs,
        CollisionNarrowphaseCollide,
        CollisionMakeContactPatch,
        QocoAssembly,
        QocoSetup,
        ClarabelAssembly,
        ClarabelSetup,
        OutputWrite,
        DynamicsAssembler_RefreshState,
        Total,
        Other
    };

    // Stats now distinguish between inclusive (wall time of the scope itself) and
    // exclusive (time excluding any child scopes).
    struct Stats {
        double totalSec{0.0};      // EXCLUSIVE time
        double inclusiveSec{0.0};  // INCLUSIVE time
        std::uint64_t calls{0};
    };

    class Scope {
       public:
        Scope(TimingManager& mgr, TimerId id) : m_mgr(&mgr), m_active(true) { m_mgr->begin(id); }
        // non-copyable
        Scope(const Scope&) = delete;
        Scope& operator=(const Scope&) = delete;
        // movable
        Scope(Scope&& other) noexcept : m_mgr(other.m_mgr), m_active(other.m_active) { other.m_active = false; }
        Scope& operator=(Scope&& other) noexcept {
            if (this != &other) {
                // end current if still active
                if (m_active && m_mgr) m_mgr->end();
                m_mgr = other.m_mgr;
                m_active = other.m_active;
                other.m_active = false;
            }
            return *this;
        }
        ~Scope() {
            if (m_active && m_mgr) m_mgr->end();
        }

       private:
        TimingManager* m_mgr{nullptr};
        bool m_active{false};
    };

    void reset() { m_stats.clear(); }

    Scope scope(TimerId id) { return Scope(*this, id); }
    Scope scope(TimerId id) const { return Scope(const_cast<TimingManager&>(*this), id); }

    const Stats& stats(TimerId id) const {
        auto it = m_stats.find(id);
        if (it != m_stats.end()) return it->second;
        static Stats empty;
        return empty;
    }

    void printBreakdown(std::ostream& os) const {
        double grandInclusive = -1.0;
        if (auto it = m_stats.find(TimerId::Total); it != m_stats.end()) grandInclusive = it->second.inclusiveSec;

        double knownExclusive = 0.0;
        for (auto& kv : m_stats) {
            if (kv.first == TimerId::Total || kv.first == TimerId::Other) continue;
            knownExclusive += kv.second.totalSec;
        }
        double grand = (grandInclusive >= 0.0) ? grandInclusive : knownExclusive;  // fallback
        if (grand <= 0.0) {
            os << "\n[Timings] No recorded timings.\n";
            return;
        }
        double otherResidual = std::max(0.0, grand - knownExclusive);

        struct Row {
            std::string name;
            double sec;
            std::uint64_t calls;
            double pct;
            double avg;
        };
        std::vector<Row> rows;
        rows.reserve(m_stats.size() + 1);
        for (auto& kv : m_stats) {
            if (kv.first == TimerId::Total || kv.first == TimerId::Other) continue;
            const auto& st = kv.second;
            const double pct = (st.totalSec / grand) * 100.0;
            const double avg = (st.calls > 0) ? (st.totalSec / (double)st.calls) : 0.0;
            rows.push_back({nameOf(kv.first), st.totalSec, st.calls, pct, avg});
        }
        if (otherResidual > 0.0 || m_stats.find(TimerId::Other) != m_stats.end()) {
            const auto calls = (m_stats.find(TimerId::Other) != m_stats.end()) ? m_stats.at(TimerId::Other).calls : 0ULL;
            const double avg = (calls > 0) ? (otherResidual / (double)calls) : 0.0;
            rows.push_back({nameOf(TimerId::Other), otherResidual, calls, (otherResidual / grand) * 100.0, avg});
        }
        std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b) { return a.sec > b.sec; });
        os << "\n[Timing Breakdown - Exclusive]\n";
        os << std::left << std::setw(24) << "Section" << std::right << std::setw(10) << "Calls" << std::setw(12) << "Excl(s)" << std::setw(10) << "%" << std::setw(12) << "Avg(µs)" << "\n";
        os << std::string(24 + 10 + 12 + 10 + 12, '-') << "\n";
        for (auto& r : rows) {
            os << std::left << std::setw(24) << r.name << std::right << std::setw(10) << r.calls << std::setw(12) << std::fixed << std::setprecision(6) << r.sec << std::setw(10)
               << std::setprecision(2) << r.pct << std::setw(12) << std::setprecision(2) << r.avg * 1e6 << "\n";
        }
        os << std::string(24 + 10 + 12 + 10 + 12, '-') << "\n";
        os << std::left << std::setw(24) << "Total (inclusive)" << std::right << std::setw(10) << "" << std::setw(12) << std::fixed << std::setprecision(4) << grand << std::setw(10) << 100.0
           << std::setw(12) << "" << "\n";

        // Hierarchy (indentation) output
        // os << "\n[Timing Hierarchy]\n";
        // double denom = grand;
        // for (auto& rid : m_roots) {
        //     auto it = m_flatNodes.find(rid);
        //     if (it != m_flatNodes.end()) printNode(os, *it->second, 0, denom);
        // }
    }

   private:
    friend class Scope;
    using clock = std::chrono::steady_clock;

    struct HierNode {
        TimerId id;
        double inclusiveSec{0.0};
        double exclusiveSec{0.0};
        std::uint64_t calls{0};
        std::unordered_map<TimerId, HierNode*> children;  // non-owning pointers to child nodes
    };

    struct Active {
        TimerId id;
        clock::time_point start;
        double childAccum{0.0};
    };

    void begin(TimerId id) { m_stack.push_back({id, clock::now(), 0.0}); }
    void end() {
        auto endTime = clock::now();
        if (m_stack.empty()) return;  // safety
        Active active = m_stack.back();
        m_stack.pop_back();
        double duration = std::chrono::duration<double>(endTime - active.start).count();
        double exclusive = duration - active.childAccum;
        auto& st = m_stats[active.id];
        st.inclusiveSec += duration;
        st.totalSec += exclusive;  // exclusive stored in totalSec (backward compat)
        ++st.calls;
        // Update hierarchy
        HierNode* node = getOrCreateNode(active.id);
        node->inclusiveSec += duration;
        node->exclusiveSec += exclusive;
        node->calls++;
        if (!m_stack.empty()) {
            // propagate inclusive duration upward
            m_stack.back().childAccum += duration;
            HierNode* parent = getOrCreateNode(m_stack.back().id);
            parent->children[active.id] = node;
        } else {
            // root node registration
            m_roots.insert(active.id);
        }
    }

    HierNode* getOrCreateNode(TimerId id) const {
        auto it = m_flatNodes.find(id);
        if (it != m_flatNodes.end()) return it->second.get();
        auto n = std::make_unique<HierNode>(HierNode{id});
        auto ptr = n.get();
        m_flatNodes[id] = std::move(n);
        return ptr;
    }

    void printNode(std::ostream& os, const HierNode& node, int depth, double grand) const {
        const double pct = (grand > 0.0) ? (node.exclusiveSec / grand) * 100.0 : 0.0;
        os << std::string(depth * 2, ' ') << "- " << nameOf(node.id) << " (incl=" << std::fixed << std::setprecision(6) << node.inclusiveSec << "s, excl=" << node.exclusiveSec << "s, "
           << std::setprecision(2) << pct << "%)\n";
        for (auto& kv : node.children) {
            if (kv.second) printNode(os, *kv.second, depth + 1, grand);
        }
    }
    static std::string nameOf(TimerId id) {
        switch (id) {
            case TimerId::CollisionBroadphase:
                return "Collision Broadphase";
            case TimerId::CollisionMatching:
                return "Collision Matching";
            case TimerId::CollisionNarrowphase:
                return "Collision Narrowphase";
            case TimerId::OutputWrite:
                return "Output Write";
            case TimerId::UpdateEntities:
                return "Update Entities";
            case TimerId::BuildAndFactorS:
                return "Build & Factor S";
            case TimerId::RebuildConstraintJacobians:
                return "Rebuild Constraints";
            case TimerId::RebuildContactJacobians:
                return "Rebuild Contacts";
            case TimerId::CollisionNarrowphaseCollide:
                return "Narrowphase Collide";
            case TimerId::Integration:
                return "Integration";
            case TimerId::CollisionMakeContactPatch:
                return "Collision Make Patches";
            case TimerId::CollisionMakeContact:
                return "Collision Make Contact";
            case TimerId::ProjectedJacobi:
                return "Projected Jacobi";
            case TimerId::ProjectedJacobiSweep:
                return "Projected Jacobi Sweep";
            case TimerId::ProjectedNewton:
                return "Projected Newton";
            case TimerId::ProjectedNewtonSetup:
                return "Projected Newton Setup";
            case TimerId::QocoSolve:
                return "QOCO Solve";
            case TimerId::ClarabelSolve:
                return "Clarabel Solve";
            case TimerId::QocoAssembly:
                return "QOCO Assembly";
            case TimerId::QocoSetup:
                return "QOCO Setup";
            case TimerId::ClarabelAssembly:
                return "Clarabel Assembly";
            case TimerId::ClarabelSetup:
                return "Clarabel Setup";
            case TimerId::DisableCollisionPairs:
                return "Disable Collision Pairs";
            case TimerId::DynamicsAssembler_RefreshState:
                return "DA Refresh State";
            case TimerId::Total:
                return "Total (all)";
            case TimerId::Other:
                return "Other";
        }
        return "Unknown";
    }
    std::unordered_map<TimerId, Stats> m_stats;
    mutable std::unordered_map<TimerId, std::unique_ptr<HierNode>> m_flatNodes;  // aggregated unique nodes per TimerId
    mutable std::unordered_set<TimerId> m_roots;                                 // root timer ids
    std::vector<Active> m_stack;                                                 // active scope stack
};

}  // namespace misc
}  // namespace cardillo