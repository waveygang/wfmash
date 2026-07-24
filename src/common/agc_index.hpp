#pragma once

/*
 * Minimal wrapper around AGC's (Assembled Genomes Compressor) public C++ API
 * (CAGCFile in agc-api.h) so wfmash can read sequences directly from a .agc
 * archive with the same operations it uses on an htslib faidx:
 *   - enumerate sequence names,
 *   - query a sequence length by name,
 *   - fetch a [start,end] (0-based, inclusive) range by name.
 *
 * The whole thing is compiled in only when the build links libagc
 * (WFMASH_HAVE_AGC, set by CMake when the AGC tree is found). When it is not
 * defined this header degrades to just is_agc_file(), and the callers'
 * #ifdef WFMASH_HAVE_AGC branches are stripped, so wfmash still builds and
 * behaves exactly as before on FASTA/FASTQ input.
 */

#include <string>

namespace agcidx {

// True if the path ends in ".agc". Always available (no AGC dependency).
inline bool is_agc_file(const std::string& filename) {
    const std::string suffix = ".agc";
    return filename.size() >= suffix.size()
        && filename.compare(filename.size() - suffix.size(), suffix.size(), suffix) == 0;
}

} // namespace agcidx

#ifdef WFMASH_HAVE_AGC

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <unordered_map>
#include <utility>
#include <vector>

#include "agc-api.h"

namespace agcidx {

/*
 * Reads one AGC archive. AGC's native identity is (sample, contig); wfmash only
 * ever knows a flat sequence name (the PAF query/target id, which for a PanSN
 * archive equals the contig name). We build an ordered (sample, contig) list and
 * a contig->sample map at open time so every flat name resolves to its sample.
 *
 * One AgcIndex owns one CAGCFile and is NOT shared between threads (the alignment
 * stage keeps one per thread, exactly like the per-thread faidx_t it replaces).
 */
class AgcIndex {
public:
    AgcIndex() = default;

    AgcIndex(const AgcIndex&) = delete;
    AgcIndex& operator=(const AgcIndex&) = delete;

    // prefetch=false keeps memory low (segments decompressed on demand), which is
    // what we want for both one-shot enumeration and per-thread random access.
    bool open(const std::string& filename, bool prefetch = false) {
        if (!agc.Open(filename, prefetch)) {
            return false;
        }
        opened = true;

        std::vector<std::string> samples;
        agc.ListSample(samples);
        for (const auto& sample : samples) {
            std::vector<std::string> contigs;
            agc.ListCtg(sample, contigs);
            for (const auto& contig : contigs) {
                ordered.emplace_back(sample, contig);
                // First occurrence of a contig name wins (mirrors impg). For a
                // PanSN archive every contig name is unique.
                contig_to_sample.emplace(contig, sample);
            }
        }
        return true;
    }

    bool is_open() const { return opened; }

    // (sample, contig) pairs in the archive's stored order.
    const std::vector<std::pair<std::string, std::string>>& sample_contigs() const {
        return ordered;
    }

    // Length of a contig by flat name; <0 if not found.
    int64_t length(const std::string& name) const {
        return agc.GetCtgLen(sample_of(name), name);
    }

    /*
     * Fetch bases [start, end], 0-based and INCLUSIVE, matching
     * faidx_fetch_seq64. end is clamped to the last base and start floored at 0 so
     * boundary requests behave identically to htslib (which silently clamps).
     * len_out receives the number of bases returned.
     */
    char* fetch_malloc(const std::string& name, int64_t start, int64_t end, int64_t& len_out) const {
        const std::string& sample = sample_of(name);
        const int64_t ctg_len = agc.GetCtgLen(sample, name);
        if (start < 0) start = 0;
        if (ctg_len > 0 && end > ctg_len - 1) end = ctg_len - 1;

        std::string buffer;
        if (end >= start) {
            agc.GetCtgSeq(sample, name, static_cast<int>(start), static_cast<int>(end), buffer);
        }

        len_out = static_cast<int64_t>(buffer.size());
        // Match htslib: caller null-terminates at len_out and frees with free().
        char* out = static_cast<char*>(std::malloc(len_out + 1));
        if (out == nullptr) {
            std::cerr << "[wfmash::agc] ERROR: out of memory fetching " << name << std::endl;
            std::exit(1);
        }
        if (len_out > 0) {
            std::memcpy(out, buffer.data(), static_cast<size_t>(len_out));
        }
        out[len_out] = '\0';
        return out;
    }

    // (contig name, length) for every contig, in archive order, using AGC metadata
    // only (GetCtgLen reads segment descriptors; it does NOT decompress the bases).
    // For the length-only code paths that otherwise decompress the whole archive.
    std::vector<std::pair<std::string, int64_t>> names_and_lengths() const {
        std::vector<std::pair<std::string, int64_t>> out;
        out.reserve(ordered.size());
        for (const auto& sc : ordered) {
            out.emplace_back(sc.second, agc.GetCtgLen(sc.first, sc.second));
        }
        return out;
    }

    // Whole contig as a std::string (used by the sequence-enumeration callback).
    std::string fetch_string(const std::string& name) const {
        const std::string& sample = sample_of(name);
        const int64_t ctg_len = agc.GetCtgLen(sample, name);
        std::string buffer;
        if (ctg_len > 0) {
            agc.GetCtgSeq(sample, name, 0, static_cast<int>(ctg_len - 1), buffer);
        }
        return buffer;
    }

private:
    const std::string& sample_of(const std::string& name) const {
        auto it = contig_to_sample.find(name);
        return it != contig_to_sample.end() ? it->second : empty_sample;
    }

    CAGCFile agc;
    bool opened = false;
    std::vector<std::pair<std::string, std::string>> ordered;
    std::unordered_map<std::string, std::string> contig_to_sample;
    std::string empty_sample; // AGC resolves an empty sample by unique contig name
};

} // namespace agcidx

#endif // WFMASH_HAVE_AGC
