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
 * Reads one AGC archive. AGC's native identity is (sample, contig), while wfmash
 * only ever knows a flat sequence name (the PAF query/target id). The mapping
 * between the two must be one-to-one, so every record gets an external name:
 *
 *   - contig name unique across all samples -> the contig name itself (the usual
 *     case for a PanSN archive, where names are already sample#hap#contig);
 *   - contig name shared by several samples -> "contig@sample", which is AGC's own
 *     query syntax, so the name still round-trips through `agc getctg`.
 *
 * Without this, an archive built the canonical AGC way (`agc create ref.fa in1.fa
 * in2.fa`, one sample per file, plain chr names) would expose several records under
 * the same name and every one of them would resolve to whichever sample was listed
 * first, silently returning the wrong sequence.
 *
 * One AgcIndex owns one CAGCFile and is NOT shared between threads (the alignment
 * stage keeps one per thread, exactly like the per-thread faidx_t it replaces).
 */
class AgcIndex {
public:
    // One archive record: the flat name wfmash sees, plus its AGC coordinates.
    struct Record {
        std::string name;
        std::string sample;
        std::string contig;
    };

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

        // Pass 1: collect (sample, contig) in archive order and count each short name.
        std::vector<std::pair<std::string, std::string>> raw;
        std::unordered_map<std::string, size_t> short_name_count;
        for (const auto& sample : samples) {
            std::vector<std::string> contigs;
            agc.ListCtg(sample, contigs);
            for (const auto& contig : contigs) {
                raw.emplace_back(sample, contig);
                ++short_name_count[short_name(contig)];
            }
        }

        // Pass 2: assign a unique external name to every record.
        ordered.reserve(raw.size());
        for (const auto& sc : raw) {
            const std::string& sample = sc.first;
            const std::string& contig = sc.second;
            const std::string base = short_name(contig);
            std::string name = short_name_count[base] > 1 ? base + "@" + sample : base;
            if (!name_to_index.emplace(name, ordered.size()).second) {
                // Only reachable if the archive itself holds a contig literally named
                // "x@y" that collides with a disambiguated name; bail out rather than
                // serve the wrong sequence.
                std::cerr << "[wfmash::agc] ERROR: duplicate sequence name '" << name
                          << "' in " << filename << std::endl;
                return false;
            }
            ordered.push_back(Record{std::move(name), sample, contig});
        }
        return true;
    }

    bool is_open() const { return opened; }

    // Every record, in the archive's stored order.
    const std::vector<Record>& records() const { return ordered; }

    // Length of a sequence by its external name; <0 if unknown.
    int64_t length(const std::string& name) const {
        const Record* rec = find(name);
        return rec ? agc.GetCtgLen(rec->sample, rec->contig) : -1;
    }

    /*
     * Fetch bases [start, end], 0-based and INCLUSIVE, matching
     * faidx_fetch_seq64. end is clamped to the last base and start floored at 0 so
     * boundary requests behave identically to htslib (which silently clamps).
     * len_out receives the number of bases returned.
     */
    char* fetch_malloc(const std::string& name, int64_t start, int64_t end, int64_t& len_out) const {
        const Record* rec = find(name);
        if (rec == nullptr) {
            std::cerr << "[wfmash::agc] ERROR: unknown sequence '" << name << "'" << std::endl;
            std::exit(1);
        }
        const int64_t ctg_len = agc.GetCtgLen(rec->sample, rec->contig);
        if (start < 0) start = 0;
        if (ctg_len > 0 && end > ctg_len - 1) end = ctg_len - 1;

        std::string buffer;
        if (end >= start) {
            agc.GetCtgSeq(rec->sample, rec->contig, static_cast<int>(start), static_cast<int>(end),
                          buffer);
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

    // (external name, length) for every record, in archive order, using AGC metadata
    // only (GetCtgLen reads segment descriptors; it does NOT decompress the bases).
    // For the length-only code paths that otherwise decompress the whole archive.
    std::vector<std::pair<std::string, int64_t>> names_and_lengths() const {
        std::vector<std::pair<std::string, int64_t>> out;
        out.reserve(ordered.size());
        for (const auto& rec : ordered) {
            out.emplace_back(rec.name, agc.GetCtgLen(rec.sample, rec.contig));
        }
        return out;
    }

    // Whole sequence as a std::string (used by the sequence-enumeration callback).
    std::string fetch_string(const std::string& name) const {
        const Record* rec = find(name);
        std::string buffer;
        if (rec == nullptr) {
            std::cerr << "[wfmash::agc] ERROR: unknown sequence '" << name << "'" << std::endl;
            return buffer;
        }
        const int64_t ctg_len = agc.GetCtgLen(rec->sample, rec->contig);
        if (ctg_len > 0) {
            agc.GetCtgSeq(rec->sample, rec->contig, 0, static_cast<int>(ctg_len - 1), buffer);
        }
        return buffer;
    }

private:
    // AGC stores the whole FASTA header line as the contig name, so ">chr1 some description"
    // comes back with its description attached. wfmash names a sequence by the first
    // whitespace-delimited token (as htslib and the FASTA reader do), and its own PAF parser
    // splits rows on whitespace, so the description must be dropped here too.
    static std::string short_name(const std::string& contig) {
        const size_t end = contig.find_first_of(" \t");
        return end == std::string::npos ? contig : contig.substr(0, end);
    }

    const Record* find(const std::string& name) const {
        auto it = name_to_index.find(name);
        return it != name_to_index.end() ? &ordered[it->second] : nullptr;
    }

    CAGCFile agc;
    bool opened = false;
    std::vector<Record> ordered;
    std::unordered_map<std::string, size_t> name_to_index;
};

} // namespace agcidx

#endif // WFMASH_HAVE_AGC
