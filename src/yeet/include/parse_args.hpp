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
    //bool align_input_paf = false;
};

void parse_args(int argc,
                char** argv,
                skch::Parameters& map_parameters,
                align::Parameters& align_parameters,
                yeet::Parameters& yeet_parameters) {

    args::ArgumentParser parser("wfmash: base-accurate alignments using mashmap2 and the wavefront algorithm");
    args::HelpFlag help(parser, "help", "display this help menu", {'h', "help"});
    args::ValueFlag<int> thread_count(parser, "N", "use this many threads during parallel steps", {'t', "threads"});
    args::Positional<std::string> target_sequence_file(parser, "target", "alignment target or reference sequence file");
    args::ValueFlag<std::string> target_sequence_file_list(parser, "targets", "alignment target file list", {'L', "target-file-list"});
    args::PositionalList<std::string> query_sequence_files(parser, "queries", "query sequences");
    args::ValueFlag<std::string> query_sequence_file_list(parser, "queries", "alignment query file list", {'Q', "query-file-list"});
    // mashmap arguments
    args::ValueFlag<int> segment_length(parser, "N", "segment length for mapping [default: 5000]", {'s', "segment-length"});
    args::ValueFlag<int> block_length_min(parser, "N", "keep mappings with at least this block length [default: 3*segment-length]", {'l', "block-length-min"});
    args::ValueFlag<int> kmer_size(parser, "N", "kmer size <= 16 [default: 16]", {'k', "kmer"});
    args::Flag no_split(parser, "no-split", "disable splitting of input sequences during mapping [enabled by default]", {'N',"no-split"});
    args::ValueFlag<float> map_pct_identity(parser, "%", "use this percent identity in the mashmap step [default: 95]", {'p', "map-pct-id"});
    args::Flag keep_low_map_pct_identity(parser, "K", "keep mappings with estimated identity below --map-pct-id=%", {'K', "keep-low-map-id"});
    args::Flag keep_low_align_pct_identity(parser, "A", "keep alignments with gap-compressed identity below --map-pct-id=%", {'O', "keep-low-align-id"});
    args::Flag no_filter(parser, "MODE", "disable mapping filtering", {'f', "no-filter"});
    args::ValueFlag<uint16_t> num_mappings_for_segments(parser, "N", "number of mappings to retain for each segment [default: 1]", {'n', "num-mappings-for-segment"});
    args::ValueFlag<uint16_t> num_mappings_for_short_seq(parser, "N", "number of mappings to retain for each sequence shorter than segment length [default: 1]", {'S', "num-mappings-for-short-seq"});
    args::Flag skip_self(parser, "", "skip self mappings when the query and target name is the same (for all-vs-all mode)", {'X', "skip-self"});
    args::ValueFlag<char> skip_prefix(parser, "C", "skip mappings when the query and target have the same prefix before the given character C", {'Y', "skip-prefix"});
    args::Flag approx_mapping(parser, "approx-map", "skip base-level alignment, producing an approximate mapping in PAF", {'m',"approx-map"});
    args::Flag no_merge(parser, "no-merge", "don't merge consecutive segment-level mappings", {'M', "no-merge"});

    args::ValueFlag<int> window_size(parser, "N", "window size for sketching. If 0, it computes the best window size applying 0 as p-value cutoff [default: automatically computed applying 1e-120 as p-value cutoff]", {'w', "window-size"});

    //args::ValueFlag<std::string> path_high_frequency_kmers(parser, "FILE", " input file containing list of high frequency kmers", {'H', "high-freq-kmers"});

    args::ValueFlag<std::string> spaced_seed_params(parser, "spaced-seed", "Params to generate spaced seeds <weight_of_seed> <number_of_seeds> <similarity> <region_length> e.g \"10 5 0.75 20\"", {'e', "spaced-seed"});

    // align parameters
    args::ValueFlag<std::string> align_input_paf(parser, "FILE", "derive precise alignments for this input PAF", {'i', "input-paf"});
    args::ValueFlag<uint16_t> wflambda_segment_length(parser, "N", "wflambda segment length: size (in bp) of segment mapped in hierarchical WFA problem [default: 256]", {'W', "wflamda-segment"});
    args::ValueFlag<int> wflambda_min_wavefront_length(parser, "N", "minimum wavefront length (width) to trigger reduction [default: 100]", {'A', "wflamda-min"});
    args::ValueFlag<int> wflambda_max_distance_threshold(parser, "N", "maximum distance that a wavefront may be behind the best wavefront [default: 100000]", {'D', "wflambda-diff"});

    //Unsupported
    //args::Flag exact_wflambda(parser, "N", "compute the exact wflambda, don't use adaptive wavefront reduction", {'xxx', "exact-wflambda"});

    // patching parameter
    args::ValueFlag<uint64_t> wflign_max_len_major(parser, "N", "maximum length to patch in the major axis [default: 512*segment-length]", {'C', "max-patch-major"});
    args::ValueFlag<uint64_t> wflign_max_len_minor(parser, "N", "maximum length to patch in the minor axis [default: 128*segment-length]", {'F', "max-patch-minor"});
    args::ValueFlag<uint16_t> wflign_erode_k(parser, "N", "maximum length of match/mismatch islands to erode before patching [default: 13]", {'E', "erode-math-mismatch"});

    // format parameters
    args::Flag emit_md_tag(parser, "N", "output the MD tag", {'d', "md-tag"});
    args::Flag sam_format(parser, "N", "output in the SAM format (PAF by default)", {'a', "sam-format"});

    // general parameters
    args::ValueFlag<std::string> tmp_base(parser, "PATH", "base name for temporary files [default: `pwd`]", {'B', "tmp-base"});
    args::Flag keep_temp_files(parser, "", "keep intermediate files generated during mapping and alignment", {'T', "keep-temp"});

    //args::Flag show_progress(parser, "show-progress", "write alignment progress to stderr", {'P', "show-progress"});
    //args::Flag verbose_debug(parser, "verbose-debug", "enable verbose debugging", {'V', "verbose-debug"});

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

    if (no_filter) {
        map_parameters.filterMode = skch::filter::NONE;
    } else {
        if (skip_self || skip_prefix) {
            // before we set skch::filter::ONETOONE here
            // but this does not provide a clear benefit in all-to-all
            // as it sometimes introduces cases of over-filtering
            map_parameters.filterMode = skch::filter::MAP;
        } else {
            map_parameters.filterMode = skch::filter::MAP;
        }
    }

    align_parameters.emit_md_tag = args::get(emit_md_tag);
    align_parameters.sam_format = args::get(sam_format);
    map_parameters.split = !args::get(no_split);

    if (align_parameters.sam_format && map_parameters.split) {
        std::cerr << "[wfmash] ERROR, skch::parseandSave, Disable splitting of input sequences (with -N) in order to enable the SAM format" << std::endl;
        exit(1);
    }

    map_parameters.mergeMappings = !args::get(no_merge);

    if (segment_length) {
        map_parameters.segLength = args::get(segment_length);
        if (map_parameters.segLength < 100) {
            std::cerr << "[wfmash] ERROR, skch::parseandSave, minimum segment length is required to be >= 100 bp." << std::endl
                      << "[wfmash] This is because Mashmap is not designed for computing short local alignments." << std::endl;
            exit(1);
        }
    } else {
        map_parameters.segLength = 5000;
    }

    if (block_length_min) {
        map_parameters.block_length_min = args::get(block_length_min);
    } else {
        map_parameters.block_length_min = 3 * map_parameters.segLength;
    }

    if (map_pct_identity) {
        if (args::get(map_pct_identity) < 70) {
            std::cerr << "[wfmash] ERROR, skch::parseandSave, minimum nucleotide identity requirement should be >= 70\%" << std::endl;
            exit(1);
        }
        map_parameters.percentageIdentity = (float) (args::get(map_pct_identity)/100.0); // scale to [0,1]
    } else {
        map_parameters.percentageIdentity = 0.95;
    }

    if (keep_low_map_pct_identity) {
        map_parameters.keep_low_pct_id = true;
    } else {
        map_parameters.keep_low_pct_id = false;
    }

    if (kmer_size) {
        map_parameters.kmerSize = args::get(kmer_size);
    } else {
        // Smaller values of k are more sensitive for divergent genomes, but lose specificity for large
        // genomes due to chance k-mer collisions. However, too large of a k-mer will reduce sensitivity
        // and so choosing the smallest k that avoids chance collisions is recommended.

        map_parameters.kmerSize = map_parameters.percentageIdentity >= 0.97 ? 19 : (map_parameters.percentageIdentity >= 0.9 ? 17 : 15);
    }

    if (spaced_seed_params) {
        auto split = [](const string& s, const string& delimiter) {
          size_t pos_start = 0, pos_end, delim_len = delimiter.length();
          string token;
          vector<string> res;

          while ((pos_end = s.find (delimiter, pos_start)) != string::npos) {
              token = s.substr (pos_start, pos_end - pos_start);
              pos_start = pos_end + delim_len;
              res.push_back (token);
          }

          res.push_back (s.substr (pos_start));
          return res;
        };

        std::string foobar = args::get(spaced_seed_params);

        // delimeters can be full colon (:) or a space
        char delimeter;
        if (foobar.find(' ') !=  std::string::npos) {
            delimeter = ' ';
        } else if (foobar.find(':') !=  std::string::npos) {
            delimeter = ':';
        } else {
            std::cerr << "[mashz] ERROR, skch::parseandSave, mashz expects either space or : for to seperate spaced seed params" << std::endl;
            exit(1);
        }

        std::string delimeter_str(1, delimeter);
        std::vector<std::string> p = split(foobar, delimeter_str);
        if (p.size() != 4) {
            std::cerr << "[mashz] ERROR, skch::parseandSave, there should be four arguments for spaced seeds" << std::endl;
            exit(1);
        }

        uint32_t seed_weight   = stoi(p[0]);
        uint32_t seed_count    = stoi(p[1]);
        float similarity       = stof(p[2]);
        uint32_t region_length = stoi(p[3]);

        // Generate an ALeS params struct
        map_parameters.use_spaced_seeds = true;
        map_parameters.spaced_seed_params = skch::ales_params{seed_weight, seed_count, similarity, region_length};
        map_parameters.kmerSize = (int) seed_weight;
    } else {
        map_parameters.use_spaced_seeds = false;
    }

    align_parameters.kmerSize = map_parameters.kmerSize;


//    if (path_high_frequency_kmers && !args::get(path_high_frequency_kmers).empty()) {
//        std::ifstream high_freq_kmers (args::get(path_high_frequency_kmers));
//
//        std::string kmer;
//        uint64_t freq;
//        while(high_freq_kmers >> kmer >> freq) {
//            if (kmer.length() != map_parameters.kmerSize) {
//                std::cerr << "[wfmash] ERROR, skch::parseandSave, high frequency k-mers length and kmerSize parameter are inconsistent." << std::endl;
//                exit(1);
//            }
//
//            map_parameters.high_freq_kmers.insert(kmer);
//        }
//        std::cerr << "[wfmash] INFO, skch::parseandSave, read " << map_parameters.high_freq_kmers.size() << " high frequency kmers." << std::endl;
//    }


    if (keep_low_align_pct_identity) {
        align_parameters.min_identity = 0; // now unused
    } else {
        align_parameters.min_identity = map_parameters.percentageIdentity; // in [0,1]
    }

    if (wflambda_segment_length) {
        align_parameters.wflambda_segment_length = args::get(wflambda_segment_length);
    } else {
        align_parameters.wflambda_segment_length = 256;
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

    if (wflign_max_len_major) {
        align_parameters.wflign_max_len_major = args::get(wflign_max_len_major);
    } else {
        align_parameters.wflign_max_len_major = map_parameters.segLength * 512;
    }

    if (wflign_max_len_minor) {
        align_parameters.wflign_max_len_minor = args::get(wflign_max_len_minor);
    } else {
        align_parameters.wflign_max_len_minor = map_parameters.segLength * 128;
    }

    if (wflign_erode_k) {
        align_parameters.wflign_erode_k = args::get(wflign_erode_k);
    } else {
        align_parameters.wflign_erode_k = map_parameters.percentageIdentity >= 0.97 ? 21 : (map_parameters.percentageIdentity >= 0.9 ? 17 : 13);
    }

    // Unsupproted
    //if (exact_wflambda) {
    //    // set exact computation of wflambda
    //    align_parameters.wflambda_min_wavefront_length = 0;
    //    align_parameters.wflambda_max_distance_threshold = 0;
    //}

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
    {
        const int ws = window_size && args::get(window_size) >= 0 ? args::get(window_size) : -1;
        if (ws > 0) {
            map_parameters.windowSize = ws;
        } else {
            // If the input window size is 0, compute the best window size using 0 as p-value cutoff
            const int windowSize = skch::Stat::recommendedWindowSize(
                    ws == 0 ? 0.0 : skch::fixed::pval_cutoff,
                    skch::fixed::confidence_interval,
                    map_parameters.kmerSize,
                    map_parameters.alphabetSize,
                    map_parameters.percentageIdentity,
                    map_parameters.segLength,
                    map_parameters.referenceSize);

            map_parameters.windowSize = std::min(
                    map_parameters.percentageIdentity >= 0.97 ? 256 : (map_parameters.percentageIdentity >= 0.9 ? 224 : 192),
                    windowSize);
        }
    }

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

    if (num_mappings_for_segments) {
        if (args::get(num_mappings_for_segments) > 0) {
            map_parameters.numMappingsForSegment = args::get(num_mappings_for_segments) ;
        } else {
            std::cerr << "[wfmash] ERROR, skch::parseandSave, the number of mappings to retain for each segment has to be grater than 0." << std::endl;
            exit(1);
        }
    } else {
        map_parameters.numMappingsForSegment = 1;
    }

    if (num_mappings_for_short_seq) {
        if (args::get(num_mappings_for_short_seq) > 0) {
            map_parameters.numMappingsForShortSequence = args::get(num_mappings_for_short_seq);
        } else {
            std::cerr << "[wfmash] ERROR, skch::parseandSave, the number of mappings to retain for each sequence shorter than segment length has to be grater than 0." << std::endl;
            exit(1);
        }
    } else {
        map_parameters.numMappingsForShortSequence = 1;
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
