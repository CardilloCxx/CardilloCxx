#pragma once

#include <string>
#include <vector>
#include "misc/types.hpp"

namespace cardillo::io {

// Accumulates (step, time, file) entries for one VTK output stream and serializes them to a
// ParaView Collection (.pvd) file, so ParaView keeps real simulation-time knowledge across the
// per-step .vtp series instead of only being able to group files by filename pattern.
class PvdWriter {
   public:
    void addEntry(int step, real_t time, const std::string& fileName);
    void write(const std::string& path) const;

   private:
    struct Entry {
        int step;
        real_t time;
        std::string fileName;
    };
    std::vector<Entry> m_entries;
};

}  // namespace cardillo::io
