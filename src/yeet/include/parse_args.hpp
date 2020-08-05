#pragma once

#include "common/args.hxx"

#include "map/include/map_parameters.hpp"
#include "map/include/map_stats.hpp"
#include "map/include/commonFunc.hpp"
#include "map/include/parseCmdArgs.hpp"

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
    args::ValueFlag<int> map_secondaries(parser, "N", "number of secondary mappings to retain in 'map' filter mode [default: 0]", {'n', "n-secondary"});
    args::Flag approx_mapping(parser, "approx-map", "skip base-level alignment, producing an approximate mapping in PAF", {'m',"approx-map"});
    // align parameters
    args::ValueFlag<std::string> align_input_paf(parser, "FILE", "derive precise alignments for this input PAF", {'i', "input-paf"});
    args::ValueFlag<float> align_pct_identity(parser, "%", "use this percent identity in the edlib step, if different than mashmap step [default: -p]", {'a', "align-pct-id"});
    args::ValueFlag<int> align_bandwidth(parser, "N", "maximum bandwidth for edlib alignment [default: 0 / computed from -p]", {'b', "align-bandwidth"});
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
            std::cerr << "ERROR, skch::parseandSave, Invalid option given for filter_mode" << std::endl;
            exit(1);
        }
    } else {
        map_parameters.filterMode = skch::filter::MAP;
    }

    map_parameters.split = args::get(no_split);
    
    if (kmer_size) {
        map_parameters.kmerSize = args::get(kmer_size);
    } else {
        map_parameters.kmerSize = 16;
    }

    if (segment_length) {
        map_parameters.segLength = args::get(segment_length);
        if (map_parameters.segLength < 500) {
            std::cerr << "ERROR, skch::parseandSave, minimum segment length is required to be >= 500 bp." << std::endl
                      << "This is because Mashmap is not designed for computing short local alignments.\n" << std::endl;
            exit(1);
        }
    } else {
        map_parameters.segLength = 5000;
    }

    if (map_pct_identity) {
        map_parameters.percentageIdentity = args::get(map_pct_identity);
        if (map_parameters.percentageIdentity < 70) {
            std::cerr << "ERROR, skch::parseandSave, minimum nucleotide identity requirement should be >= 70\%\n" << std::endl;
            exit(1);
        }
    } else {
        map_parameters.percentageIdentity = 85;
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
        if (tmp_base) {
            map_parameters.outFileName = temp_file::create(args::get(tmp_base));
        } else {
            map_parameters.outFileName = temp_file::create();
        }
        align_parameters.mashmapPafFile = map_parameters.outFileName;
        align_parameters.pafOutputFile = "/dev/stdout";
    }

    if (map_secondaries) {
        map_parameters.secondaryToKeep = args::get(map_secondaries);
    } else {
        map_parameters.secondaryToKeep = 0;
    }

    skch::printCmdOptions(map_parameters);

    //Check if files are valid
    skch::validateInputFiles(map_parameters.querySequences, map_parameters.refSequences);

    // set up alignment parameters

    

}

}
