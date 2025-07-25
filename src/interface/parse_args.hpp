#pragma once

#include <unistd.h>
#include <limits.h>
#include <regex>

#include "common/args.hxx"

#include "map/include/map_parameters.hpp"
#include "map/include/map_stats.hpp"
#include "map/include/commonFunc.hpp"

#include "align/include/align_parameters.hpp"

#include "interface/temp_file.hpp"
#include "common/utils.hpp"

#include "wfmash_git_version.hpp"

// If the WFMASH_GIT_VERSION doesn't exist at all, define a placeholder
#ifndef WFMASH_GIT_VERSION
#define WFMASH_GIT_VERSION "not-from-git"
#endif

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

    args::ArgumentParser parser("");
    parser.helpParams.width = 100;
    parser.helpParams.showTerminator = false;
    parser.helpParams.progindent = 0;
    parser.helpParams.descriptionindent = 0;
    parser.helpParams.flagindent = 2;
    parser.helpParams.helpindent = 35;
    parser.helpParams.eachgroupindent = 2;
    parser.helpParams.optionsString = "";

    args::Group options_group(parser, "");

    // SEQUENCES
    args::Positional<std::string> target_sequence_file(options_group, "target.fa", "target sequences (required)");
    args::Positional<std::string> query_sequence_file(options_group, "query.fa", "query sequences (default: self-map)");
    // INDEXING
    args::Group indexing_opts(options_group, "INDEXING:");
    args::ValueFlag<std::string> write_index(indexing_opts, "FILE", "save index to FILE", {'W', "write-index"});
    args::ValueFlag<std::string> read_index(indexing_opts, "FILE", "load index from FILE", {'I', "read-index"});
    args::ValueFlag<std::string> index_by(indexing_opts, "SIZE", "target batch size for indexing [4G]", {'b', "batch"});

    // MINMERS
    args::Group minmer_opts(options_group, "MINMERS:");
    args::ValueFlag<int> kmer_size(minmer_opts, "INT", "k-mer size [15]", {'k', "kmer-size"});
    args::ValueFlag<int64_t> sketch_size(minmer_opts, "INT", "sketch size [auto]", {'s', "sketch-size"});
    args::ValueFlag<std::string> window_size(minmer_opts, "INT", "window size [1k]", {'w', "window-size"});

    // MAPPING
    args::Group mapping_opts(options_group, "MAPPING:");
    args::Flag approx_mapping(mapping_opts, "", "output mappings only (no alignment)", {'m', "approx-mapping"});
    args::ValueFlag<std::string> map_pct_identity(mapping_opts, "FLOAT|aniXX[+/-N]", "minimum identity % or ANI preset (default: ani25-5)", {'p', "map-pct-id"});
    args::ValueFlag<int> ani_sketch_size(mapping_opts, "INT", "sketch size for ANI estimation [1000]", {"ani-sketch-size"});
    args::ValueFlag<uint32_t> num_mappings(mapping_opts, "INT", "mappings per segment [1]", {'n', "mappings"});
    args::ValueFlag<std::string> block_length(mapping_opts, "INT", "minimum block length [0]", {'l', "block-length"});
    args::ValueFlag<std::string> chain_jump(mapping_opts, "INT", "chain jump (gap) [2k]", {'c', "chain-jump"});
    args::ValueFlag<std::string> max_mapping_length(mapping_opts, "INT", "max mapping length [50k]", {'P', "max-length"});
    args::Flag no_split(mapping_opts, "", "map each sequence in one piece", {'N', "no-split"});

    // FILTERING
    args::Group filtering_opts(options_group, "FILTERING:");
    args::Flag no_filter(filtering_opts, "", "disable all filtering", {'f', "no-filter"});
    args::Flag no_merge(filtering_opts, "", "keep all fragment mappings", {'M', "no-merge"});
    args::Flag one_to_one(filtering_opts, "", "best mapping per query AND target", {'o', "one-to-one"});
    args::ValueFlag<double> overlap_threshold(filtering_opts, "FLOAT", "max overlap ratio [1.0]", {'O', "overlap"});
    args::ValueFlag<double> map_sparsification(filtering_opts, "FLOAT", "keep this fraction of mappings [1.0]", {'x', "sparsify"});
    args::ValueFlag<std::string> hg_filter(filtering_opts, "n,Δ,conf", "hypergeometric filter [1.0,0.0,99.9]", {"hg-filter"});
    args::ValueFlag<int> min_hits(filtering_opts, "INT", "min hits for L1 filtering [3]", {'H', "l1-hits"});
    args::ValueFlag<double> max_kmer_freq(filtering_opts, "FLOAT", "filter high-freq minimizers [0.0002]", {'F', "filter-freq"});

    // SCAFFOLDING
    args::Group scaffold_opts(options_group, "SCAFFOLDING:");
    args::ValueFlag<int> scaffold_mass(scaffold_opts, "INT", "scaffold mass (min segments) [5]", {'S', "scaffold-mass"});
    args::ValueFlag<std::string> scaffold_dist(scaffold_opts, "INT", "max scaffold distance [100k]", {'D', "scaffold-dist"});
    args::ValueFlag<std::string> scaffold_jump(scaffold_opts, "INT", "scaffold jump (gap) [100k]", {'j', "scaffold-jump"});
    args::ValueFlag<std::string> scaffold_output(scaffold_opts, "FILE", "output scaffold mappings to FILE", {"scaffold-out"});

    // SELECTION
    args::Group selection_opts(options_group, "SELECTION:");
    args::ValueFlag<char> skip_prefix(selection_opts, "C", "group by prefix delimiter [#]", {'Y', "group-prefix"});
    args::Flag enable_self_mappings(selection_opts, "", "include self mappings", {'X', "self-maps"});
    args::Flag lower_triangular(selection_opts, "", "only map seq_i to seq_j if i>j", {'L', "lower-triangular"});
    args::ValueFlag<std::string> target_prefix(selection_opts, "pfx", "target name prefix filter", {'T', "target-prefix"});
    args::ValueFlag<std::string> target_list(selection_opts, "FILE", "file with target sequence names", {'R', "target-list"});
    args::ValueFlag<std::string> query_prefix(selection_opts, "pfxs", "query name prefix filter", {'Q', "query-prefix"});
    args::ValueFlag<std::string> query_list(selection_opts, "FILE", "file with query sequence names", {'A', "query-list"});

    // ALIGNMENT
    args::Group alignment_opts(options_group, "ALIGNMENT:");
    args::ValueFlag<std::string> input_mapping(alignment_opts, "FILE", "align from PAF file", {'i', "align-paf"});
    args::ValueFlag<std::string> target_padding(alignment_opts, "INT", "target padding [segment-length]", {'E', "target-padding"});
    args::ValueFlag<std::string> query_padding(alignment_opts, "INT", "query padding [segment-length]", {'U', "query-padding"});
    args::ValueFlag<std::string> wfa_params(alignment_opts, "m,go1,ge1,go2,ge2", "gap costs [5,8,2,24,1]", {'g', "wfa-params"});

    // OUTPUT
    args::Group output_opts(options_group, "OUTPUT:");
    args::Flag sam_format(output_opts, "", "SAM format (default: PAF)", {'a', "sam"});
    args::Flag emit_md_tag(output_opts, "", "include MD tag", {'d', "md-tag"});

    // SYSTEM
    args::Group system_opts(options_group, "SYSTEM:");
    args::ValueFlag<int> thread_count(system_opts, "INT", "threads [1]", {'t', "threads"});
    args::ValueFlag<std::string> tmp_base(system_opts, "PATH", "temp file directory [pwd]", {'B', "tmp-base"});
    args::Flag keep_temp_files(system_opts, "", "retain temporary files", {'Z', "keep-temp"});
    args::Flag quiet(system_opts, "", "disable progress output", {"quiet"});
    args::Flag streaming_minhash(system_opts, "", "use streaming MinHash algorithm (experimental)", {"streaming-minhash"});
    args::Flag version(system_opts, "", "version info", {'v', "version"});
    args::HelpFlag help(system_opts, "", "show help", {'h', "help"});

#ifdef WFA_PNG_TSV_TIMING
    args::Group debugging_opts(parser, "[ Debugging Options ]");
    args::ValueFlag<std::string> prefix_wavefront_info_in_tsv(parser, "PREFIX", " write wavefronts' information for each alignment in TSV format files with this PREFIX", {'G', "tsv"});
    args::ValueFlag<std::string> prefix_wavefront_plot_in_png(parser, "PREFIX", " write wavefronts' plot for each alignment in PNG format files with this PREFIX", {'u', "prefix-png"});
    args::ValueFlag<uint64_t> wfplot_max_size(parser, "N", "max size of the wfplot [default: 1500]", {'z', "wfplot-max-size"});
    args::ValueFlag<std::string> path_patching_info_in_tsv(parser, "FILE", " write patching information for each alignment in TSV format in FILE", {"path-patching-tsv"});
#endif


    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        exit(0);
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        exit(1);
    }

    if (version) {
        std::cerr << WFMASH_GIT_VERSION << std::endl;
        exit(0);
    }

    if (argc==1 || !target_sequence_file) {
        std::cout << parser;
        exit(1);
    }

    map_parameters.skip_self = !args::get(enable_self_mappings);
    map_parameters.lower_triangular = lower_triangular ? args::get(lower_triangular) : false;
    map_parameters.keep_low_pct_id = true;

    map_parameters.use_progress_bar = !args::get(quiet);
    align_parameters.use_progress_bar = !args::get(quiet);
    map_parameters.use_streaming_minhash = args::get(streaming_minhash);
    
    if (ani_sketch_size) {
        map_parameters.ani_sketch_size = args::get(ani_sketch_size);
    }

    if (skip_prefix) {
        map_parameters.prefix_delim = args::get(skip_prefix);
        map_parameters.skip_prefix = map_parameters.prefix_delim != '\0';
    } else {
        map_parameters.prefix_delim = '#';
        map_parameters.skip_prefix = true;
    }

	if (target_list) {
		map_parameters.target_list = args::get(target_list);
	}

	if (target_prefix) {
		map_parameters.target_prefix = args::get(target_prefix);
	}

	if (query_list) {
		map_parameters.query_list = args::get(query_list);
	}
	
	if (query_prefix) {
		map_parameters.query_prefix = skch::CommonFunc::split(args::get(query_prefix), ',');
	}
	
    map_parameters.refSequences.push_back(args::get(target_sequence_file));
    align_parameters.refSequences.push_back(args::get(target_sequence_file));

    map_parameters.referenceSize = skch::CommonFunc::getReferenceSize(map_parameters.refSequences);

    if (query_sequence_file) {
        map_parameters.querySequences.push_back(args::get(query_sequence_file));
        align_parameters.querySequences.push_back(args::get(query_sequence_file));
    }


    // If there are no queries, go in all-vs-all mode with the sequences specified in `target_sequence_file`
    if (map_parameters.querySequences.empty()) {
        std::cerr << "[wfmash] Performing all-vs-all mapping including self mappings." << std::endl;
        map_parameters.querySequences.push_back(map_parameters.refSequences.back());
        align_parameters.querySequences.push_back(align_parameters.refSequences.back());
    }

    map_parameters.alphabetSize = 4;

    if (no_filter) {
        map_parameters.filterMode = skch::filter::NONE;
    } else {
        if (one_to_one) {
            map_parameters.filterMode = skch::filter::ONETOONE;
        } else {
            map_parameters.filterMode = skch::filter::MAP;
        }
    }

    if (map_sparsification) {
        if (args::get(map_sparsification) == 1) {
            // overflows
            map_parameters.sparsity_hash_threshold = std::numeric_limits<uint64_t>::max();
        } else {
            map_parameters.sparsity_hash_threshold
                = args::get(map_sparsification) * std::numeric_limits<uint64_t>::max();
        }
    } else {
        map_parameters.sparsity_hash_threshold
            = std::numeric_limits<uint64_t>::max();
    }

    args::ValueFlag<std::string> wfa_score_params(alignment_opts, "MISMATCH,GAP,EXT", "WFA scoring parameters [2,3,1]", {"wfa-params"});
    if (!args::get(wfa_score_params).empty()) {
        const std::vector<std::string> params_str = skch::CommonFunc::split(args::get(wfa_score_params), ',');
        if (params_str.size() != 3) {
            std::cerr << "[wfmash] ERROR error: 3 scoring parameters must be given to --wfa-params"
                      << std::endl;
            exit(1);
        }

        std::vector<int> params(params_str.size());
        std::transform(params_str.begin(), params_str.end(), params.begin(),
                       [](const std::string &s) { return std::stoi(s); });

        align_parameters.wfa_mismatch_score = params[0];
        align_parameters.wfa_gap_opening_score = params[1];
        align_parameters.wfa_gap_extension_score = params[2];
    } else {
        align_parameters.wfa_mismatch_score = 2;
        align_parameters.wfa_gap_opening_score = 3;
        align_parameters.wfa_gap_extension_score = 1;
    }

    if (!args::get(wfa_params).empty()) {
        const std::vector<std::string> params_str = skch::CommonFunc::split(args::get(wfa_params), ',');
        if (params_str.size() != 5) {
            std::cerr << "[wfmash] ERROR error: 5 scoring parameters must be given to --wfa-params"
                      << std::endl;
            exit(1);
        }

        std::vector<int> params(params_str.size());
        std::transform(params_str.begin(), params_str.end(), params.begin(),
                       [](const std::string &s) { return std::stoi(s); });

        align_parameters.wfa_patching_mismatch_score = params[0];
        align_parameters.wfa_patching_gap_opening_score1 = params[1];
        align_parameters.wfa_patching_gap_extension_score1 = params[2];
        align_parameters.wfa_patching_gap_opening_score2 = params[3];
        align_parameters.wfa_patching_gap_extension_score2 = params[4];
    } else {
        align_parameters.wfa_patching_mismatch_score = 5;
        align_parameters.wfa_patching_gap_opening_score1 = 8;
        align_parameters.wfa_patching_gap_extension_score1 = 2;
        align_parameters.wfa_patching_gap_opening_score2 = 24;
        align_parameters.wfa_patching_gap_extension_score2 = 1;
    }

    align_parameters.wflign_mismatch_score = 2;
    align_parameters.wflign_gap_opening_score = 3;
    align_parameters.wflign_gap_extension_score = 1;

    align_parameters.wflign_max_mash_dist = -1;
    align_parameters.wflign_min_wavefront_length = 1024;
    align_parameters.wflign_max_distance_threshold = -1;

    align_parameters.emit_md_tag = args::get(emit_md_tag);
    align_parameters.sam_format = args::get(sam_format);
    align_parameters.no_seq_in_sam = false;
    align_parameters.disable_chain_patching = false;
    args::Flag force_wflign(alignment_opts, "", "force WFlign alignment", {"force-wflign"});
    align_parameters.force_wflign = args::get(force_wflign);
    map_parameters.split = !args::get(no_split);
    map_parameters.dropRand = false;//ToFix: !args::get(keep_ties);
    align_parameters.split = !args::get(no_split);

    map_parameters.mergeMappings = !args::get(no_merge);

    if (window_size) {
        const int64_t s = wfmash::handy_parameter(args::get(window_size));

        if (s <= 0) {
            std::cerr << "[wfmash] ERROR, skch::parseandSave, window size has to be a float value greater than 0." << std::endl;
            exit(1);
        }

        if (s < 100) {
            std::cerr << "[wfmash] ERROR, skch::parseandSave, minimum window size is required to be >= 100 bp." << std::endl
                      << "[wfmash] This is because Mashmap is not designed for computing short local alignments." << std::endl;
            exit(1);
        }

        if (!approx_mapping && s > 10000) {
            std::cerr << "[wfmash] ERROR: window size (-w) must be <= 10kb when running alignment." << std::endl
                      << "[wfmash] For larger values, use -m/--approx-mapping to generate mappings," << std::endl
                      << "[wfmash] then align them with: wfmash ... -i mappings.paf" << std::endl;
            exit(1);
        }
        map_parameters.windowLength = s;
    } else {
        map_parameters.windowLength = 1000; // Default 1k
    }

    if (map_pct_identity) {
        std::string pct_id_str = args::get(map_pct_identity);
        
        // Check if it's an ANI preset (e.g., ani25, ani50-10, ani75+5)
        std::regex ani_pattern("^ani(\\d+)([+-]\\d+)?$");
        std::smatch matches;
        
        if (std::regex_match(pct_id_str, matches, ani_pattern)) {
            map_parameters.auto_pct_identity = true;
            map_parameters.ani_percentile = std::stoi(matches[1]);
            
            // Validate percentile
            if (map_parameters.ani_percentile < 1 || map_parameters.ani_percentile > 99) {
                std::cerr << "[wfmash] ERROR: ANI percentile must be between 1 and 99, got: " 
                          << map_parameters.ani_percentile << std::endl;
                exit(1);
            }
            
            // Parse adjustment if present
            if (matches[2].matched) {
                map_parameters.ani_adjustment = std::stof(matches[2].str());
            } else {
                map_parameters.ani_adjustment = 0.0;
            }
            
            map_parameters.percentageIdentity = skch::fixed::percentage_identity; // Will be overridden
        } else if (pct_id_str == "auto") {
            // Legacy support for "auto" = ani25
            map_parameters.auto_pct_identity = true;
            map_parameters.ani_percentile = 25;
            map_parameters.ani_adjustment = 0.0;
            map_parameters.percentageIdentity = skch::fixed::percentage_identity;
        } else {
            // Try to parse as a float percentage
            try {
                float pct_id = std::stof(pct_id_str);
                if (pct_id < 50) {
                    std::cerr << "[wfmash] ERROR: minimum nucleotide identity requirement should be >= 50%." << std::endl;
                    exit(1);
                }
                map_parameters.percentageIdentity = pct_id / 100.0; // scale to [0,1]
                map_parameters.auto_pct_identity = false;
            } catch (const std::exception& e) {
                std::cerr << "[wfmash] ERROR: Invalid value for -p/--map-pct-id: " << pct_id_str << std::endl;
                std::cerr << "[wfmash] Expected a float value or ANI preset (e.g., ani25, ani50-10, ani75+5)" << std::endl;
                exit(1);
            }
        }
    } else {
        // No -p flag provided, use ani25-5 by default
        map_parameters.percentageIdentity = skch::fixed::percentage_identity; // Will be overridden
        // auto_pct_identity, ani_percentile=25, ani_adjustment=-5.0 already set in Parameters struct
    }

    if (block_length) {
        const int64_t l = wfmash::handy_parameter(args::get(block_length));

        if (l < 0) {
            std::cerr << "[wfmash] ERROR, skch::parseandSave, min block length has to be a float value greater than or equal to 0." << std::endl;
            exit(1);
        }

        if (!approx_mapping && l > 30000) {
            std::cerr << "[wfmash] ERROR: block length (-l) must be <= 30kb when running alignment." << std::endl
                      << "[wfmash] For larger values, use -m/--approx-mapping to generate mappings," << std::endl
                      << "[wfmash] then align them with: wfmash ... -i mappings.paf" << std::endl;
            exit(1);
        }
        map_parameters.block_length = l;
    } else {
        map_parameters.block_length = 0; // Default 0 (no filtering)
    }

    if (chain_jump) {
        const int64_t l = wfmash::handy_parameter(args::get(chain_jump));
        if (l < 0) {
            std::cerr << "[wfmash] ERROR, skch::parseandSave, chain jump has to be a float value greater than or equal to 0." << std::endl;
            exit(1);
        }
        map_parameters.chain_gap = l;
        align_parameters.chain_gap = l;
    } else {
        map_parameters.chain_gap = 2000;
        align_parameters.chain_gap = 2000;
    }

    // Parse scaffold jump (gap)
    if (scaffold_jump) {
        map_parameters.scaffold_gap = wfmash::handy_parameter(args::get(scaffold_jump));
        if (map_parameters.scaffold_gap < 0) {
            std::cerr << "[wfmash] ERROR: Invalid scaffold jump" << std::endl;
            exit(1);
        }
    } else {
        map_parameters.scaffold_gap = 100000; // 100k default
    }

    // Parse scaffold distance threshold
    if (scaffold_dist) {
        map_parameters.scaffold_max_deviation = wfmash::handy_parameter(args::get(scaffold_dist));
        if (map_parameters.scaffold_max_deviation < 0) {
            std::cerr << "[wfmash] ERROR: Invalid scaffold distance" << std::endl;
            exit(1);
        }
    } else {
        map_parameters.scaffold_max_deviation = 100000; // 100k default
    }

    // Parse scaffold mass (minimum segments)
    if (scaffold_mass) {
        int min_segs = args::get(scaffold_mass);
        if (min_segs <= 0) {
            std::cerr << "[wfmash] ERROR: Invalid scaffold mass" << std::endl;
            exit(1);
        }
        // Convert segments to length: segments * window_length
        map_parameters.scaffold_min_length = min_segs * map_parameters.windowLength;
    } else {
        // Default: 5 segments of 1kb each = 5kb
        map_parameters.scaffold_min_length = 5000;
    }

    if (scaffold_output) {
        map_parameters.scaffold_output_file = args::get(scaffold_output);
    }

    if (max_mapping_length) {
        const int64_t l = args::get(max_mapping_length) == "inf" ? std::numeric_limits<int64_t>::max()
            : wfmash::handy_parameter(args::get(max_mapping_length));
        if (l <= 0) {
            std::cerr << "[wfmash] ERROR: max mapping length must be greater than 0." << std::endl;
            exit(1);
        }
        // Removed restriction on max mapping length
        map_parameters.max_mapping_length = l;
    } else {
        map_parameters.max_mapping_length = 50000;
    }

    if (map_parameters.windowLength >= map_parameters.max_mapping_length) {
        std::cerr << "[wfmash] ERROR, skch::parseandSave, window size should not be larger than max mapping length." << std::endl;
        exit(1);
    }
    if (map_parameters.block_length >= map_parameters.max_mapping_length) {
        std::cerr << "[wfmash] ERROR, skch::parseandSave, block length should not be larger than max mapping length." << std::endl;
        exit(1);
    }

    if (overlap_threshold) {
        map_parameters.overlap_threshold = args::get(overlap_threshold);
    } else {
        map_parameters.overlap_threshold = 1.0;
    }

    if (kmer_size) {
        map_parameters.kmerSize = args::get(kmer_size);
    } else {
        // Smaller values of k are more sensitive for divergent genomes, but lose specificity for large
        // genomes due to chance k-mer collisions. However, too large of a k-mer will reduce sensitivity
        // and so choosing the smallest k that avoids chance collisions is recommended.
        /*
        map_parameters.kmerSize = (map_parameters.percentageIdentity >= 0.97 ? 18 :
                                  (map_parameters.percentageIdentity >= 0.9 ? 17 : 15));
        */
        map_parameters.kmerSize = 15;
    }

    //if (spaced_seed_params) {
        //const std::string foobar = args::get(spaced_seed_params);

        //// delimeters can be full colon (:) or a space
        //char delimeter;
        //if (foobar.find(' ') !=  std::string::npos) {
            //delimeter = ' ';
        //} else if (foobar.find(':') !=  std::string::npos) {
            //delimeter = ':';
        //} else {
            //std::cerr << "[wfmash] ERROR, skch::parseandSave, wfmash expects either space or : for to seperate spaced seed params." << std::endl;
            //exit(1);
        //}

        //const std::vector<std::string> p = skch::CommonFunc::split(foobar, delimeter);
        //if (p.size() != 4) {
            //std::cerr << "[wfmash] ERROR, skch::parseandSave, there should be four arguments for spaced seeds." << std::endl;
            //exit(1);
        //}

        //const uint32_t seed_weight   = stoi(p[0]);
        //const uint32_t seed_count    = stoi(p[1]);
        //const float similarity       = stof(p[2]);
        //const uint32_t region_length = stoi(p[3]);

        //// Generate an ALeS params struct
        //map_parameters.use_spaced_seeds = true;
        //map_parameters.spaced_seed_params = skch::ales_params{seed_weight, seed_count, similarity, region_length};
        //map_parameters.kmerSize = (int) seed_weight;
    //} else {
        //map_parameters.use_spaced_seeds = false;
    //}
    map_parameters.use_spaced_seeds = false;

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

    align_parameters.min_identity = 0; // disabled

    args::ValueFlag<int> wflambda_segment_length(alignment_opts, "N", "WFlambda segment length [256]", {"wflambda-segment"});
    if (wflambda_segment_length) {
        align_parameters.wflambda_segment_length = args::get(wflambda_segment_length);
    } else {
        align_parameters.wflambda_segment_length = 256;
    }

    align_parameters.wflign_max_len_major = map_parameters.windowLength * 512;
    align_parameters.wflign_max_len_minor = map_parameters.windowLength * 128;
    align_parameters.wflign_erode_k = -1; // will trigger estimation based on sequence divergence
    align_parameters.wflign_min_inv_patch_len = 23;
    align_parameters.wflign_max_patching_score = 0; // will trigger estimation based on gap penalties and sequence length

    if (target_padding) {
        const int64_t p = wfmash::handy_parameter(args::get(target_padding));
        if (p < 0) {
            std::cerr << "[wfmash] ERROR: target padding must be >= 0" << std::endl;
            exit(1);
        }
        align_parameters.target_padding = p;
    } else {
        // Default to half window length, capped at 5000
        align_parameters.target_padding = std::min(map_parameters.windowLength, (int64_t)5000);
    }
    
    if (query_padding) {
        const int64_t p = wfmash::handy_parameter(args::get(query_padding));
        if (p < 0) {
            std::cerr << "[wfmash] ERROR: query padding must be >= 0" << std::endl;
            exit(1);
        }
        align_parameters.query_padding = p;
    } else {
        // Default to half window length, capped at 5000
        align_parameters.query_padding = std::min(map_parameters.windowLength, (int64_t)5000);
    }

    if (thread_count) {
        map_parameters.threads = args::get(thread_count);
        align_parameters.threads = args::get(thread_count);
    } else {
        map_parameters.threads = 1;
        align_parameters.threads = 1;
    }
    // disable multi-fasta processing due to the memory inefficiency of samtools faidx readers
    // which require us to duplicate the in-memory indexes of large files for each thread
    // if aligner exhaustion is a problem, we could enable this
    align_parameters.multithread_fasta_input = false;

    // Compute optimal window size for sketching
    {
        const int64_t ss = sketch_size && args::get(sketch_size) >= 0 ? args::get(sketch_size) : -1;
        if (ss > 0) {
            map_parameters.sketchSize = ss;
        } else {
            const double md = 1 - map_parameters.percentageIdentity;
            double dens = 0.02 * (1 + (md / 0.1));
            map_parameters.sketchSize = dens * (map_parameters.windowLength - map_parameters.kmerSize);
        }
    }

    // Validate that sketch size <= window size
    if (map_parameters.sketchSize > map_parameters.windowLength) {
        std::cerr << "[wfmash] ERROR: sketch size (-s) must be less than or equal to window size (-w)." << std::endl;
        std::cerr << "[wfmash] Current values: sketch size = " << map_parameters.sketchSize 
                  << ", window size = " << map_parameters.windowLength << std::endl;
        exit(1);
    }

    map_parameters.kmerComplexityThreshold = 0;

    args::ValueFlag<double> hg_numerator(mapping_opts, "FLOAT", "hypergeometric filter numerator [1.0]", {"hg-numerator"});
    if (hg_numerator) {
        double value = args::get(hg_numerator);
        if (value < 1.0) {
            std::cerr << "[wfmash] ERROR: hg-numerator must be >= 1.0." << std::endl;
            exit(1);
        }
        map_parameters.hgNumerator = value;
    } else {
        map_parameters.hgNumerator = 1.0;  // Default value
    }

    // Create sequence ID manager for getting sequence info
    std::unique_ptr<skch::SequenceIdManager> idManager = std::make_unique<skch::SequenceIdManager>(
        map_parameters.querySequences,
        map_parameters.refSequences,
        std::vector<std::string>{map_parameters.query_prefix},
        std::vector<std::string>{map_parameters.target_prefix},
        std::string(1, map_parameters.prefix_delim),
        map_parameters.query_list,
        map_parameters.target_list);

    // Get target sequence names
    std::vector<std::string> targetSequenceNames = idManager->getTargetSequenceNames();

    // Estimate total unique k-mers using information theoretic approach
    map_parameters.estimatedUniqueKmers = skch::CommonFunc::estimateUniqueKmers(
        map_parameters.totalReferenceSize, 
        map_parameters.kmerSize
    );

    // Calculate total reference size from actual sequence lengths
    map_parameters.totalReferenceSize = 0;
    // Calculate total reference size silently
    for (const auto& seqName : targetSequenceNames) {
        skch::seqno_t seqId = idManager->getSequenceId(seqName);
        skch::offset_t seqLen = idManager->getSequenceLength(seqId);
        map_parameters.totalReferenceSize += seqLen;
    }

    map_parameters.filterLengthMismatches = true;

    // Parse hypergeometric filter parameters
    map_parameters.stage1_topANI_filter = true;
    map_parameters.stage2_full_scan = true;
    
    if (hg_filter) {
        std::string hg_params = args::get(hg_filter);
        std::vector<std::string> params = skch::CommonFunc::split(hg_params, ',');
        if (params.size() != 3) {
            std::cerr << "[wfmash] ERROR: hypergeometric filter requires 3 comma-separated values: numerator,ani-diff,confidence" << std::endl;
            exit(1);
        }
        // Parse numerator
        map_parameters.hgNumerator = std::stod(params[0]);
        if (map_parameters.hgNumerator < 1.0) {
            std::cerr << "[wfmash] ERROR: hg-filter numerator must be >= 1.0" << std::endl;
            exit(1);
        }
        // Parse ANI difference
        map_parameters.ANIDiff = std::stod(params[1]) / 100.0;
        // Parse confidence
        map_parameters.ANIDiffConf = std::stod(params[2]) / 100.0;
    } else {
        // Use defaults
        map_parameters.hgNumerator = 1.0;
        map_parameters.ANIDiff = skch::fixed::ANIDiff;
        map_parameters.ANIDiffConf = skch::fixed::ANIDiffConf;
    }

    if (min_hits) {
        map_parameters.minimum_hits = args::get(min_hits);
    } else {
        map_parameters.minimum_hits = 3; // default minimum
    }

    if (max_kmer_freq) {
        map_parameters.max_kmer_freq = args::get(max_kmer_freq);
    } else {
        map_parameters.max_kmer_freq = 0.0002; // default filter fraction
    }

    //if (window_minimizers) {
        //map_parameters.world_minimizers = false;
    //} else {
        //map_parameters.world_minimizers = true;
    //}

    if (write_index) {
        map_parameters.indexFilename = args::get(write_index);
        map_parameters.overwrite_index = true;
        map_parameters.create_index_only = true;
    } else if (read_index) {
        map_parameters.indexFilename = args::get(read_index);
        map_parameters.overwrite_index = false;
        map_parameters.create_index_only = false;
    } else {
        map_parameters.indexFilename = "";
        map_parameters.overwrite_index = false;
        map_parameters.create_index_only = false;
    }

    if (index_by) {
        const int64_t index_size = wfmash::handy_parameter(args::get(index_by));
        if (index_size < 0) {
            std::cerr << "[wfmash] ERROR, skch::parseandSave, index-by size must be a positive integer." << std::endl;
            exit(1);
        }
        map_parameters.index_by_size = static_cast<int64_t>(index_size);
    } else {
        map_parameters.index_by_size = std::numeric_limits<int64_t>::max(); // Default to indexing all sequences
    }

    if (approx_mapping) {
        map_parameters.outFileName = "/dev/stdout";
        yeet_parameters.approx_mapping = true;
    } else {
        yeet_parameters.approx_mapping = false;

        if (tmp_base) {
            temp_file::set_dir(args::get(tmp_base));
        } else {
            char cwd[PATH_MAX];
            if (getcwd(cwd, sizeof(cwd)) != NULL) {
                temp_file::set_dir(std::string(cwd));
            } else {
                // Handle error: getcwd() failed
                std::cerr << "[wfmash] ERROR, skch::parseandSave, problem in getting the current directory." << std::endl;
                exit(1);
            }
        }

        if (input_mapping) {
            // directly use the input PAF file
            yeet_parameters.remapping = true;
            map_parameters.outFileName = args::get(input_mapping);
            align_parameters.mashmapPafFile = args::get(input_mapping);
        } else {
            // make a temporary mapping file
            map_parameters.outFileName = temp_file::create();
            align_parameters.mashmapPafFile = map_parameters.outFileName;
        }
        align_parameters.pafOutputFile = "/dev/stdout";
    }

#ifdef WFA_PNG_TSV_TIMING
    align_parameters.tsvOutputPrefix = (prefix_wavefront_info_in_tsv && !args::get(prefix_wavefront_info_in_tsv).empty())
            ? args::get(prefix_wavefront_info_in_tsv)
            : "";

    align_parameters.path_patching_info_in_tsv = (path_patching_info_in_tsv && !args::get(path_patching_info_in_tsv).empty())
                                               ? args::get(path_patching_info_in_tsv)
                                               : "";

    // wfplotting
    if (prefix_wavefront_plot_in_png) {
        align_parameters.prefix_wavefront_plot_in_png = args::get(prefix_wavefront_plot_in_png);
    }
    if (wfplot_max_size) {
        align_parameters.wfplot_max_size = args::get(wfplot_max_size);
    } else {
        align_parameters.wfplot_max_size = 1500;
    }
#endif

    if (num_mappings) {
        if (args::get(num_mappings) > 0) {
            map_parameters.numMappingsForSegment = args::get(num_mappings);
        } else {
            std::cerr << "[wfmash] ERROR: the number of mappings to retain (-n) must be greater than 0." << std::endl;
            exit(1);
        }
    } else {
        map_parameters.numMappingsForSegment = 1;
    }

    map_parameters.numMappingsForShortSequence = 1;

	map_parameters.legacy_output = false;

    //Check if files are valid
    skch::validateInputFiles(map_parameters.querySequences, map_parameters.refSequences);

    std::cerr << "[wfmash] Parameters: k=" << map_parameters.kmerSize 
              << ", s=" << map_parameters.sketchSize
              << ", w=" << map_parameters.windowLength << (map_parameters.split ? " (split)" : "")
              << ", l=" << map_parameters.block_length
              << ", c=" << map_parameters.chain_gap
              << ", P=" << map_parameters.max_mapping_length
              << ", j=" << map_parameters.scaffold_gap 
              << ", D=" << map_parameters.scaffold_max_deviation
              << ", S=" << (map_parameters.scaffold_min_length / map_parameters.windowLength)
              << ", n=" << map_parameters.numMappingsForSegment
              << ", p=";
    
    // Show the ANI preset if auto_pct_identity is enabled
    if (map_parameters.auto_pct_identity) {
        std::cerr << "ani" << map_parameters.ani_percentile;
        if (map_parameters.ani_adjustment != 0) {
            std::cerr << std::showpos << map_parameters.ani_adjustment << std::noshowpos;
        }
    } else {
        std::cerr << std::fixed << std::setprecision(0) << map_parameters.percentageIdentity * 100 << "%";
    }
    
    std::cerr
              << ", t=" << map_parameters.threads
              << ", b=" << map_parameters.index_by_size << std::endl;
    std::cerr << "[wfmash] Filters: " << (map_parameters.skip_self ? "skip-self" : "no-skip-self")
              << ", hg(Δ=" << map_parameters.ANIDiff << ",conf=" << map_parameters.ANIDiffConf << ")"
              << ", mode=" << map_parameters.filterMode << " (1=map,2=1-to-1,3=none)" << std::endl;
    std::cerr << "[wfmash] Output: " << map_parameters.outFileName << std::endl;

    temp_file::set_keep_temp(args::get(keep_temp_files));
}

}
