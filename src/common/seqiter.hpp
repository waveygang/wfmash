#pragma once

#include <string>
#include <functional>
#include <cassert>
#include <unordered_set>
#include "gzstream.h"
#include <htslib/faidx.h>

namespace seqiter {

bool fai_index_exists(const std::string& filename) {
  // Check if .fai file exists
  std::ifstream f(filename + ".fai");
  return f.good();
}

void for_each_seq_in_file(
    const std::string& filename,
    const std::unordered_set<std::string>& keep_seq,
    const std::string& keep_prefix,
    const std::function<void(const std::string&, const std::string&)>& func) {

    if ((!keep_seq.empty() || !keep_prefix.empty())
          && fai_index_exists(filename)) {
        // Use index
        // Handle keep_prefix
        std::unordered_set<std::string> found_seq;
        auto faid = fai_load(filename.c_str());
        if (!keep_prefix.empty()) {
            // iterate over lines in the fai file
            std::string line;
            std::ifstream in(filename + ".fai");
            while (std::getline(in, line)) {
                std::string name = line.substr(0, line.find("\t"));
                char* seq;
                int64_t len = 0;
                if (strncmp(name.c_str(), keep_prefix.c_str(), keep_prefix.size()) == 0
                    || keep_seq.find(name) != keep_seq.end()) {
                    found_seq.insert(name);
                    seq = faidx_fetch_seq64(
                        faid, name.c_str(), 0, INT_MAX, &len);
                    func(name, std::string(seq, len));
                    free(seq);
                } else {
                    func(name, "");
                }
            }
        }
        // Handle keep_seq
        for (const auto& name : keep_seq) {
            if (found_seq.find(name) == found_seq.end())
            {
                std::cerr << "[wfmash::for_each_seq_in_file] could not fetch " << name << " from index" << std::endl;
            }
        }
        fai_destroy(faid); // Free FAI index
    } else {
        // no index available
        // detect file type
        bool input_is_fasta = false;
        bool input_is_fastq = false;
        std::string line;
        igzstream in(filename.c_str());
        std::getline(in, line);
        if (line[0] == '>') {
            input_is_fasta = true;
        } else if (line[0] == '@') {
            input_is_fastq = true;
        } else {
            std::cerr << "[wfmash::for_each_seq_in_file] unknown file format given to seqiter" << std::endl;
            assert(false);
            exit(1);
        }
        if (input_is_fasta) {
            while (in.good()) {
                std::string name = line.substr(1, line.find(" ")-1);
                std::string seq;
                bool keep = (keep_prefix.empty() || name.substr(0, keep_prefix.length()) == keep_prefix)
                    && (keep_seq.empty() || keep_seq.find(name) != keep_seq.end());
                while (std::getline(in, line)) {
                    if (line[0] == '>') {
                        // this is the header of the next sequence
                        break;
                    } else {
                        if (keep) {
                            seq.append(line);
                        }
                    }
                }
                func(name, seq);
            }
        } else if (input_is_fastq) {
            while (in.good()) {
                std::string name = line.substr(1, line.find(" ")-1);
                std::string seq;
                bool keep = (keep_prefix.empty() || name.substr(0, keep_prefix.length()) == keep_prefix)
                    && (keep_seq.empty() || keep_seq.find(name) != keep_seq.end());
                std::getline(in, seq); // sequence
                std::getline(in, line); // delimiter
                std::getline(in, line); // quality
                std::getline(in, line); // next header
                func(name, keep ? seq : "");
            }
        }
    }
}

void for_each_seq_in_faidx_t(
    faidx_t* fai,
    const std::vector<std::string>& seq_names,
    const std::function<void(const std::string&, const std::string&)>& func) {
    for (const auto& seq_name : seq_names) {
        int len;
        char* seq = fai_fetch(fai, seq_name.c_str(), &len);
        if (seq != nullptr) {
            func(seq_name, std::string(seq));
            free(seq);
        }
    }
}

void for_each_seq_in_file(
    const std::string& filename,
    const std::vector<std::string>& seq_names,
    const std::function<void(const std::string&, const std::string&)>& func) {
    faidx_t* fai = fai_load(filename.c_str());
    for_each_seq_in_faidx_t(fai, seq_names, func);
    fai_destroy(fai);
}
	
void for_each_seq_in_file_filtered(
    const std::string& filename,
    const std::vector<std::string>& query_prefix,
    const std::unordered_set<std::string>& query_list,
    const std::function<void(const std::string&, const std::string&)>& func) {
    faidx_t* fai = fai_load(filename.c_str());
    if (fai == nullptr) {
        std::cerr << "Error: Failed to load FASTA index for file " << filename << std::endl;
        return;
    }

    std::vector<std::string> query_seq_names;
    int num_seqs = faidx_nseq(fai);
    for (int i = 0; i < num_seqs; i++) {
        const char* seq_name = faidx_iseq(fai, i);
        bool keep = false;
        for (const auto& prefix : query_prefix) {
            if (strncmp(seq_name, prefix.c_str(), prefix.size()) == 0) {
                keep = true;
                break;
            }
        }
        if (query_list.empty() || query_list.count(seq_name)) {
            keep = true;
        }
        if (keep) {
            query_seq_names.push_back(seq_name);
        }
    }

    for_each_seq_in_faidx_t(
        fai,
        query_seq_names,
        func);

    fai_destroy(fai);
}

} // namespace seqiter
