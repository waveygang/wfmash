#pragma once

#include <unistd.h>

#include "common/args.hxx"

#include "map/include/map_parameters.hpp"
#include "map/include/map_stats.hpp"
#include "map/include/commonFunc.hpp"

#include "align/include/align_parameters.hpp"

#include "yeet/include/temp_file.hpp"

namespace yeet {

struct Parameters {
    bool approx_mapping = false;
    bool remapping = false;
};

void parse_args(int argc,
                char** argv,
                skch::Parameters& map_parameters,
                align::Parameters& align_parameters,
                yeet::Parameters& yeet_parameters) {

    args::ArgumentParser parser("edyeet: base-accurate alignments using edlib and mashmap2");
    args::HelpFlag help(parser, "help", "display this help menu", {'h', "help"});
    args::ValueFlag<uint64_t> thread_count(parser, "N", "use this many threads during parallel steps", {'t', "threads"});
    args::Positional<std::string> target_sequence_file(parser, "target", "alignment target or reference sequence file");
    args::ValueFlag<std::string> target_sequence_file_list(parser, "targets", "alignment target file list", {'L', "target-file-list"});
    args::PositionalList<std::string> query_sequence_files(parser, "queries", "query sequences");
    args::ValueFlag<std::string> query_sequence_file_list(parser, "queries", "alignment query file list", {'Q', "query-file-list"});
    // mashmap arguments
    args::ValueFlag<uint64_t> segment_length(parser, "N", "segment length for mapping [default: 5000]", {'s', "segment-length"});
    args::Flag no_split(parser, "no-split", "disable splitting of input sequences during mapping [enabled by default]", {'N',"no-split"});
    args::ValueFlag<float> map_pct_identity(parser, "%", "use this percent identity in the mashmap step [default: 85]", {'p', "map-pct-id"});
    args::ValueFlag<int> kmer_size(parser, "N", "kmer size <= 16 [default: 16]", {'k', "kmer"});
    args::ValueFlag<std::string> map_filter_mode(parser, "MODE", "filter mode for map step, either 'map', 'one-to-one', or 'none' [default: map]", {'f', "map-filter-mode"});
    args::ValueFlag<int> map_secondaries(parser, "N", "number of secondary mappings to retain in 'map' filter mode (total number of mappings is this + 1) [default: 0]", {'n', "n-secondary"});
    args::ValueFlag<int> map_short_secondaries(parser, "N", "number of secondary mappings to retain for sequences shorter than segment length [default: 0]", {'S', "n-short-secondary"});
    args::Flag skip_self(parser, "", "skip self mappings when the query and target name is the same (for all-vs-all mode)", {'X', "skip-self"});
    args::ValueFlag<char> skip_prefix(parser, "C", "skip mappings when the query and target have the same prefix before the given character C", {'Y', "skip-prefix"});
    args::Flag approx_mapping(parser, "approx-map", "skip base-level alignment, producing an approximate mapping in PAF", {'m',"approx-map"});
    args::Flag merge_mappings(parser, "merge-map", "merge consecutive segment-level mappings (can slow alignment phase)", {'M', "merge-mappings"});
    // align parameters
    args::ValueFlag<std::string> align_input_paf(parser, "FILE", "derive precise alignments for this input PAF", {'i', "input-paf"});
    args::ValueFlag<float> align_pct_identity(parser, "%", "use this percent identity in the edlib step, if different than mashmap step [default: -p]", {'a', "align-pct-id"});
    args::ValueFlag<int> wflambda_segment_length(parser, "N", "wflambda segment length: size (in bp) of segment mapped in hierarchical WFA problem [default: 1000]", {'W', "wflamda-segment"});
    args::ValueFlag<int> wflambda_min_wavefront_length(parser, "N", "minimum wavefront length (width) to trigger reduction [default: 100]", {'A', "wflamda-min"});
    args::ValueFlag<int> wflambda_max_distance_threshold(parser, "N", "maximum distance that a wavefront may be behind the best wavefront [default: 100000]", {'D', "wflambda-diff"});
    args::Flag exact_wflambda(parser, "N", "compute the exact wflambda, don't use adaptive wavefront reduction", {'E', "exact-wflambda"});

    // general parameters
    args::ValueFlag<std::string> tmp_base(parser, "PATH", "base name for temporary files [default: `pwd`]", {'B', "tmp-base"});
    args::Flag keep_temp_files(parser, "", "keep intermediate files generated during mapping and alignment", {'T', "keep-temp"});
    args::Flag show_progress(parser, "show-progress", "write alignment progress to stderr", {'P', "show-progress"});
    args::Flag verbose_debug(parser, "verbose-debug", "enable verbose debugging", {'V', "verbose-debug"});

    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        exit(0);
        //return; // 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        exit(1);
        //return; // 1;
    }
    if (argc==1) {
        std::cout << parser;
        exit(1);
        //return; // 1;
    }


    if (target_sequence_file) {
        map_parameters.refSequences.push_back(args::get(target_sequence_file));
        align_parameters.refSequences.push_back(args::get(target_sequence_file));
    }
    if (target_sequence_file_list) {
        skch::parseFileList(args::get(target_sequence_file_list), map_parameters.refSequences);
        skch::parseFileList(args::get(target_sequence_file_list), align_parameters.refSequences);
    }
    map_parameters.referenceSize = skch::CommonFunc::getReferenceSize(map_parameters.refSequences);

    if (query_sequence_files) {
        for (auto& q : args::get(query_sequence_files)) {
            map_parameters.querySequences.push_back(q);
            align_parameters.querySequences.push_back(q);
        }
    }
    if (query_sequence_file_list) {
        skch::parseFileList(args::get(query_sequence_file_list), map_parameters.querySequences);
        skch::parseFileList(args::get(query_sequence_file_list), align_parameters.querySequences);
    }
    
    map_parameters.alphabetSize = 4;

    if (map_filter_mode) {
        auto& filter_input = args::get(map_filter_mode);
        if (filter_input == "map") map_parameters.filterMode = skch::filter::MAP;
        else if (filter_input == "one-to-one") map_parameters.filterMode = skch::filter::ONETOONE;
        else if (filter_input == "none") map_parameters.filterMode = skch::filter::NONE;
        else 
        {
            std::cerr << "[edyeet] ERROR, skch::parseandSave, Invalid option given for filter_mode" << std::endl;
            exit(1);
        }
    } else {
        map_parameters.filterMode = skch::filter::MAP;
    }

    map_parameters.split = !args::get(no_split);
    map_parameters.mergeMappings = args::get(merge_mappings);
    
    if (kmer_size) {
        map_parameters.kmerSize = args::get(kmer_size);
    } else {
        map_parameters.kmerSize = 16;
    }

    if (segment_length) {
        map_parameters.segLength = args::get(segment_length);
        if (map_parameters.segLength < 200) {
            std::cerr << "[edyeet] ERROR, skch::parseandSave, minimum segment length is required to be >= 200 bp." << std::endl
                      << "[edyeet] This is because Mashmap is not designed for computing short local alignments." << std::endl;
            exit(1);
        }
    } else {
        map_parameters.segLength = 5000;
    }

    if (map_pct_identity) {
        map_parameters.percentageIdentity = args::get(map_pct_identity);
        if (map_parameters.percentageIdentity < 70) {
            std::cerr << "[edyeet] ERROR, skch::parseandSave, minimum nucleotide identity requirement should be >= 70\%" << std::endl;
            exit(1);
        }
    } else {
        map_parameters.percentageIdentity = 85;
    }

    if (align_pct_identity) {
        align_parameters.percentageIdentity = args::get(align_pct_identity);
    } else {
        align_parameters.percentageIdentity = map_parameters.percentageIdentity;
    }

    if (wflambda_segment_length) {
        align_parameters.wflambda_segment_length = args::get(wflambda_segment_length);
    } else {
        align_parameters.wflambda_segment_length = 1000;
    }

    if (wflambda_min_wavefront_length) {
        align_parameters.wflambda_min_wavefront_length = args::get(wflambda_min_wavefront_length);
    } else {
        align_parameters.wflambda_min_wavefront_length = 100;
    }

    if (wflambda_max_distance_threshold) {
        align_parameters.wflambda_max_distance_threshold = args::get(wflambda_max_distance_threshold);
    } else {
        align_parameters.wflambda_max_distance_threshold = 100000;
    }
    align_parameters.wflambda_max_distance_threshold /= (align_parameters.wflambda_segment_length / 2); // set relative to WFA matrix

    if (exact_wflambda) {
        // set exact computation of wflambda
        align_parameters.wflambda_min_wavefront_length = 0;
        align_parameters.wflambda_max_distance_threshold = 0;
    }

    if (thread_count) {
        map_parameters.threads = args::get(thread_count);
        align_parameters.threads = args::get(thread_count);
    } else {
        map_parameters.threads = 1;
        align_parameters.threads = 1;
    }

    /*
     * Compute window size for sketching
     */

    //Compute optimal window size
    map_parameters.windowSize = skch::Stat::recommendedWindowSize(skch::fixed::pval_cutoff,
                                                                  map_parameters.kmerSize,
                                                                  map_parameters.alphabetSize,
                                                                  map_parameters.percentageIdentity,
                                                                  map_parameters.segLength,
                                                                  map_parameters.referenceSize);


    if (approx_mapping) {
        map_parameters.outFileName = "/dev/stdout";
        yeet_parameters.approx_mapping = true;
    } else {
        yeet_parameters.approx_mapping = false;
        if (align_input_paf) {
            yeet_parameters.remapping = true;
            align_parameters.mashmapPafFile = args::get(align_input_paf);
        } else {
            if (tmp_base) {
                temp_file::set_dir(args::get(tmp_base));
                map_parameters.outFileName = temp_file::create();
            } else {
                char* cwd = get_current_dir_name();
                temp_file::set_dir(std::string(cwd));
                free(cwd);
                map_parameters.outFileName = temp_file::create();
            }
            align_parameters.mashmapPafFile = map_parameters.outFileName;
        }
        align_parameters.pafOutputFile = "/dev/stdout";
    }

    if (map_secondaries) {
        map_parameters.secondaryToKeep = args::get(map_secondaries);
    } else {
        map_parameters.secondaryToKeep = 0;
    }

    if (map_short_secondaries) {
        map_parameters.shortSecondaryToKeep = args::get(map_short_secondaries);
    } else {
        map_parameters.shortSecondaryToKeep = 0;
    }

    if (skip_self) {
        map_parameters.skip_self = true;
    } else {
        map_parameters.skip_self = false;
    }

    if (skip_prefix) {
        map_parameters.skip_prefix = true;
        map_parameters.prefix_delim = args::get(skip_prefix);
    } else {
        map_parameters.skip_prefix = false;
        map_parameters.prefix_delim = '\0';
    }

    //Check if files are valid
    skch::validateInputFiles(map_parameters.querySequences, map_parameters.refSequences);

    temp_file::set_keep_temp(args::get(keep_temp_files));

}

}
