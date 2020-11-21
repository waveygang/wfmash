#include <iostream>
#include "args.hxx"
#include "seqiter.hpp"
#include "wflign.hpp"

void parse_file_list(const std::string& file,
                     std::vector<std::string>& files) {
    std::ifstream in(file);
    if (in.fail()) {
        std::cerr << "[wflign::main] error, could not open " << file << std::endl;
        exit(1);
    }
    std::string line;
    while (std::getline(in, line)) {
        files.push_back(line);
    }
}

int main(int argc, char** argv) {

    args::ArgumentParser parser("wflign: we-flyin wavefront-guided aligner");
    args::HelpFlag help(parser, "help", "display this help menu", {'h', "help"});
    //args::ValueFlag<uint64_t> thread_count(parser, "N", "use this many threads during parallel steps", {'t', "threads"});
    args::Positional<std::string> target_sequence_file(parser, "target", "alignment target or reference sequence file");
    args::ValueFlag<std::string> target_sequence_file_list(parser, "targets", "alignment target file list", {'L', "target-file-list"});
    args::PositionalList<std::string> query_sequence_files(parser, "queries", "query sequences");
    args::ValueFlag<std::string> query_sequence_file_list(parser, "queries", "alignment query file list", {'Q', "query-file-list"});
    args::ValueFlag<uint64_t> p_segment_length(parser, "N", "segment length for aligning [default: 1000]", {'s', "segment-length"});
    args::ValueFlag<float> min_pct_identity(parser, "%", "only emit alignments above this percent identity [default: 0]", {'I', "min-pct-id"});
    args::ValueFlag<int> wf_min(parser, "N", "WF_min: minimum length of a wavefront to trigger reduction [default: 100]", {'l', "wf-min"});
    args::ValueFlag<int> wf_diff(parser, "N", "WF_diff: maximum distance that a wavefront may be behind the best wavefront to not be reduced [default: 50]", {'d', "wf-diff"});
    args::Flag exact_wfa(parser, "N", "compute the exact WFA, don't use adaptive wavefront reduction", {'e', "exact-wfa"});

    // general parameters
    //args::Flag show_progress(parser, "show-progress", "write alignment progress to stderr", {'P', "show-progress"});
    //args::Flag verbose_debug(parser, "verbose-debug", "enable verbose debugging", {'V', "verbose-debug"});

    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    if (argc==1) {
        std::cout << parser;
        return 1;
    }

    if (!target_sequence_file && !target_sequence_file_list || !query_sequence_files && !query_sequence_file_list) {
        std::cerr << "[wflign::main] error: both a target and query must be provided" << std::endl;
        return 1;
    }

    std::vector<std::string> targets;
    if (target_sequence_file) {
        targets.push_back(args::get(target_sequence_file));
    }
    if (target_sequence_file_list) {
        parse_file_list(args::get(target_sequence_file_list), targets);
    }
    std::vector<std::string> queries;
    for (auto& q : args::get(query_sequence_files)) {
        queries.push_back(q);
    }
    if (query_sequence_file_list) {
        parse_file_list(args::get(query_sequence_file_list), queries);
    }

    uint64_t segment_length = p_segment_length ? args::get(p_segment_length) : 1000;
    uint64_t min_wavefront_length = wf_min ? args::get(wf_min) : 100;
    uint64_t max_distance_threshold = wf_diff ? args::get(wf_diff) : 200;
    float min_identity = min_pct_identity ? args::get(min_pct_identity) / 100 : 0;

    // exact WFA is triggered by setting the reduction parameters to 0
    if (args::get(exact_wfa)) {
        min_wavefront_length = 0;
        max_distance_threshold = 0;
    }

    // simple all-vs-all mapping for testing
    for (auto& query_file : queries) {
        seqiter::for_each_seq_in_file(
            query_file,
            [&](const std::string& qname,
                const std::string& qseq) {
                for (auto& target_file : targets) {
                    seqiter::for_each_seq_in_file(
                        target_file,
                        [&](const std::string& tname,
                            const std::string& tseq) {
                            wflign::wflign_affine_wavefront(
                                std::cout,
                                qname, qseq.c_str(), qseq.size(), 0, qseq.size(),
                                false, // query is rev
                                tname, tseq.c_str(), tseq.size(), 0, tseq.size(),
                                segment_length,
                                min_identity,
                                min_wavefront_length,
                                max_distance_threshold);
                        });
                }
            });
    }

    return 0;
}
