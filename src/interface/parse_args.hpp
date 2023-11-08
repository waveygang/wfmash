#pragma once

#include <unistd.h>

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
    //bool align_input_paf = false;
};

int64_t handy_parameter(const std::string& value) {
    auto is_a_float = [](const std::string s) {
        return !s.empty() && s.find_first_not_of("0123456789.") == std::string::npos && std::count(s.begin(), s.end(), '.') < 2;
    };

    uint64_t str_len = value.length();
    uint8_t exp = 0;
    if (value[str_len-1] == 'k' || value[str_len-1] == 'K') {
        exp = 3;
        --str_len;
    } else if (value[str_len-1] == 'm' || value[str_len-1] == 'M') {
        exp = 6;
        --str_len;
    } else if (value[str_len-1] == 'g' || value[str_len-1] == 'G') {
        exp = 9;
        --str_len;
    }

    const std::string tmp = value.substr(0, str_len);
    return is_a_float(tmp) ? (int)(stof(tmp) * pow(10, exp)) : -1;
}

void parse_args(int argc,
                char** argv,
                skch::Parameters& map_parameters,
                align::Parameters& align_parameters,
                yeet::Parameters& yeet_parameters) {

    args::ArgumentParser parser("wfmash: a pangenome-scale aligner, " + std::string(WFMASH_GIT_VERSION));
    parser.helpParams.width = 100;
    parser.helpParams.showTerminator = false;

    args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
    args::Positional<std::string> target_sequence_file(mandatory_opts, "target", "alignment target/reference sequence file");

    args::Group io_opts(parser, "[ Files IO Options ]");
    args::PositionalList<std::string> query_sequence_files(io_opts, "queries", "query sequences file");
    args::ValueFlag<std::string> query_sequence_file_list(io_opts, "queries", "alignment queries files list", {'Q', "query-file-list"});

    args::Group mapping_opts(parser, "[ Mapping Options ]");
    args::ValueFlag<float> map_pct_identity(mapping_opts, "%", "percent identity in the mashmap step [default: 90]", {'p', "map-pct-id"});
    args::ValueFlag<std::string> segment_length(mapping_opts, "N", "segment seed length for mapping [default: 5k]", {'s', "segment-length"});
    args::ValueFlag<std::string> block_length(mapping_opts, "N", "keep merged mappings supported by homologies of this total length [default: 5*segment-length]", {'l', "block-length"});
    args::ValueFlag<uint32_t> num_mappings_for_segments(mapping_opts, "N", "number of mappings to retain for each segment [default: 1]", {'n', "num-mappings-for-segment"});
    args::ValueFlag<uint32_t> num_mappings_for_short_seq(mapping_opts, "N", "number of mappings to retain for each sequence shorter than segment length [default: 1]", {'S', "num-mappings-for-short-seq"});
    args::ValueFlag<int> kmer_size(mapping_opts, "N", "kmer size [default: 19]", {'k', "kmer"});
    args::ValueFlag<float> kmer_pct_threshold(mapping_opts, "%", "ignore the top % most-frequent kmers [default: 0.001]", {'H', "kmer-threshold"});
	args::Flag lower_triangular(mapping_opts, "", "only map shorter sequences against longer", {'L', "lower-triangular"});
    args::Flag skip_self(mapping_opts, "", "skip self mappings when the query and target name is the same (for all-vs-all mode)", {'X', "skip-self"});
    args::ValueFlag<char> skip_prefix(mapping_opts, "C", "skip mappings when the query and target have the same prefix before the last occurrence of the given character C", {'Y', "skip-prefix"});
	args::ValueFlag<std::string> target_prefix(mapping_opts, "pfx", "use only targets whose name starts with this prefix", {'P', "target-prefix"});
	args::ValueFlag<std::string> target_list(mapping_opts, "FILE", "file containing list of target sequence names to use", {'A', "target-list"});
    args::Flag approx_mapping(mapping_opts, "approx-map", "skip base-level alignment, producing an approximate mapping in PAF", {'m',"approx-map"});
    args::Flag no_split(mapping_opts, "no-split", "disable splitting of input sequences during mapping [default: enabled]", {'N',"no-split"});
    args::ValueFlag<std::string> chain_gap(mapping_opts, "N", "chain mappings closer than this distance in query and target, sets approximate maximum variant length detectable in alignment [default: 20k]", {'c', "chain-gap"});
    args::Flag drop_low_map_pct_identity(mapping_opts, "K", "drop mappings with estimated identity below --map-pct-id=%", {'K', "drop-low-map-id"});
    args::Flag no_filter(mapping_opts, "MODE", "disable mapping filtering", {'f', "no-filter"});
    args::ValueFlag<double> map_sparsification(mapping_opts, "FACTOR", "keep this fraction of mappings", {'x', "sparsify-mappings"});
    args::Flag keep_same_score(mapping_opts, "", "keep all mappings with equal score even if it results in more than n mappings", {'D', "keep-same-score"});
    args::ValueFlag<int64_t> sketch_size(mapping_opts, "N", "sketch size for sketching.", {'w', "sketch-size"});
    args::ValueFlag<double> kmer_complexity(mapping_opts, "F", "Drop segments w/ predicted kmer complexity below this cutoff. Kmer complexity defined as #kmers / (s - k + 1)", {'J', "kmer-complexity"});
    args::Flag no_hg_filter(mapping_opts, "", "Don't use the hypergeometric filtering and instead use the MashMap2 first pass filtering.", {'1', "no-hg-filter"});
    args::ValueFlag<double> hg_filter_ani_diff(mapping_opts, "%", "Filter out mappings unlikely to be this ANI less than the best mapping [default: 0.0]", {'2', "hg-filter-ani-diff"});
    args::ValueFlag<double> hg_filter_conf(mapping_opts, "%", "Confidence value for the hypergeometric filtering [default: 99.9%]", {'3', "hg-filter-conf"});
    //args::Flag window_minimizers(mapping_opts, "", "Use window minimizers rather than world minimizers", {'U', "window-minimizers"});
    //args::ValueFlag<std::string> path_high_frequency_kmers(mapping_opts, "FILE", " input file containing list of high frequency kmers", {'H', "high-freq-kmers"});
    //args::ValueFlag<std::string> spaced_seed_params(mapping_opts, "spaced-seeds", "Params to generate spaced seeds <weight_of_seed> <number_of_seeds> <similarity> <region_length> e.g \"10 5 0.75 20\"", {'e', "spaced-seeds"});
    args::Flag no_merge(mapping_opts, "no-merge", "don't merge consecutive segment-level mappings", {'M', "no-merge"});

    args::Group alignment_opts(parser, "[ Alignment Options ]");
    args::ValueFlag<std::string> align_input_paf(alignment_opts, "FILE", "derive precise alignments for this input PAF", {'i', "input-paf"});
    args::Flag invert_filtering(alignment_opts, "A", "if an input PAF is specified, remove alignments with gap-compressed identity below --map-pct-id x 0.8, else keep all alignments "
                                                   "[default: if an input PAF is specified, keep all alignments, else remove alignments with gap-compressed identity below --map-pct-id x 0.8]",
                                {'O', "invert-filtering"});
    args::ValueFlag<uint16_t> wflambda_segment_length(alignment_opts, "N", "wflambda segment length: size (in bp) of segment mapped in hierarchical WFA problem [default: 256]", {'W', "wflamda-segment"});
    args::ValueFlag<std::string> wfa_score_params(alignment_opts, "mismatch,gap1,ext1",
                                            "score parameters for the wfa alignment (affine); match score is fixed at 0 [default: 4,6,1]",//, if 4 values then gaps are affine, if 6 values then gaps are convex",
                                            {'g', "wfa-params"});
    //wflign parameters
    args::ValueFlag<std::string> wflign_score_params(alignment_opts, "mismatch,gap1,ext1",
                                                       "score parameters for the wflign alignment (affine); match score is fixed at 0 [default: 4,6,1]",//, if 4 then gaps are affine, if 6 then gaps are convex [default: 1,4,6,2,26,1]",
                                                       {'G', "wflign-params"});
    args::ValueFlag<float> wflign_max_mash_dist(alignment_opts, "N", "maximum mash distance to perform the alignment in a wflambda segment [default: adaptive with respect to the estimated identity]", {'b', "max-mash-dist"});
    args::ValueFlag<int> wflign_min_wavefront_length(alignment_opts, "N", "min wavefront length for heuristic WFlign [default: 1024]", {'j', "wflign-min-wf-len"});
    args::ValueFlag<int> wflign_max_distance_threshold(alignment_opts, "N", "max distance threshold for heuristic WFlign [default: 2048/(estimated_identity^2)]", {'q', "wflign-max-distance"});

    // patching parameter
    args::ValueFlag<std::string> wflign_max_len_major(alignment_opts, "N", "maximum length to patch in the major axis [default: 512*segment-length]", {'C', "max-patch-major"});
    args::ValueFlag<std::string> wflign_max_len_minor(alignment_opts, "N", "maximum length to patch in the minor axis [default: 128*segment-length]", {'F', "max-patch-minor"});
    args::ValueFlag<int> wflign_erode_k(alignment_opts, "N", "maximum length of match/mismatch islands to erode before patching [default: adaptive]", {'E', "erode-match-mismatch"});
    args::ValueFlag<int> wflign_max_patching_score(alignment_opts, "N", "maximum score allowed when patching [default: adaptive with respect to gap penalties and sequence length]", {"max-patching-score"});

    args::Group output_opts(parser, "[ Output Format Options ]");
    // format parameters
    args::Flag emit_md_tag(output_opts, "N", "output the MD tag", {'d', "md-tag"});
    // sam format
    args::Flag sam_format(output_opts, "N", "output in the SAM format (PAF by default)", {'a', "sam-format"});
    args::Flag no_seq_in_sam(output_opts, "N", "do not fill the sequence field in the SAM format", {'q', "no-seq-in-sam"});

    args::Group general_opts(parser, "[ General Options ]");
    args::ValueFlag<std::string> tmp_base(general_opts, "PATH", "base name for temporary files [default: `pwd`]", {'B', "tmp-base"});
    args::Flag keep_temp_files(general_opts, "", "keep intermediate files", {'Z', "keep-temp"});
    //args::Flag show_progress(general_opts, "show-progress", "write alignment progress to stderr", {'P', "show-progress"});

#ifdef WFA_PNG_AND_TSV
    args::Group debugging_opts(parser, "[ Debugging Options ]");
    args::ValueFlag<std::string> prefix_wavefront_info_in_tsv(parser, "PREFIX", " write wavefronts' information for each alignment in TSV format files with this PREFIX", {'T', "tsv"});
    args::ValueFlag<std::string> prefix_wavefront_plot_in_png(parser, "PREFIX", " write wavefronts' plot for each alignment in PNG format files with this PREFIX", {'u', "prefix-png"});
    args::ValueFlag<uint64_t> wfplot_max_size(parser, "N", "max size of the wfplot [default: 1500]", {'z', "wfplot-max-size"});
    args::ValueFlag<std::string> path_patching_info_in_tsv(parser, "FILE", " write patching information for each alignment in TSV format in FILE", {"path-patching-tsv"});
#endif

    args::Group threading_opts(parser, "[ Threading ]");
    args::ValueFlag<int> thread_count(threading_opts, "N", "use this many threads during parallel steps", {'t', "threads"});

    args::Group program_info_opts(parser, "[ Program Information ]");
    args::Flag version(program_info_opts, "version", "show version number and github commit hash", {'v', "version"});
    args::HelpFlag help(program_info_opts, "help", "display this help menu", {'h', "help"});

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

    if (version) {
        std::cerr << WFMASH_GIT_VERSION << std::endl;
        exit(0);
    }

    if (skip_self) {
        map_parameters.skip_self = true;
    } else {
        map_parameters.skip_self = false;
    }

	if (lower_triangular) {
		map_parameters.lower_triangular = true;
	} else {
		map_parameters.lower_triangular = false;
	}

    if (skip_prefix) {
        map_parameters.skip_prefix = true;
        map_parameters.prefix_delim = args::get(skip_prefix);
    } else {
        map_parameters.skip_prefix = false;
        map_parameters.prefix_delim = '\0';
    }

	if (target_list) {
		map_parameters.target_list = args::get(target_list);
	}

	if (target_prefix) {
		map_parameters.target_prefix = args::get(target_prefix);
	}
	
    if (target_sequence_file) {
        map_parameters.refSequences.push_back(args::get(target_sequence_file));
        align_parameters.refSequences.push_back(args::get(target_sequence_file));
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

	if (target_sequence_file && map_parameters.querySequences.empty()
		&& map_parameters.refSequences.size() == 1
		&& !map_parameters.lower_triangular
		&& map_parameters.target_list.empty()
		&& map_parameters.target_prefix.empty()) {
		std::cerr << "[wfmash] Warning: Detected single file all-vs-all mapping with no other options. "
				  << "Consider adding -L, --lower-triangular for efficiency." << std::endl;
	}

    // If there are no queries, go in all-vs-all mode with the sequences specified in `target_sequence_file`
    if (target_sequence_file && map_parameters.querySequences.empty()) {
        map_parameters.skip_self = true;
        std::cerr << "[mashmap] Skipping self mappings for single file all-vs-all mapping." << std::endl;
        map_parameters.querySequences.push_back(map_parameters.refSequences.back());
        align_parameters.querySequences.push_back(align_parameters.refSequences.back());
    }

    map_parameters.alphabetSize = 4;

    if (no_filter) {
        map_parameters.filterMode = skch::filter::NONE;
    } else {
        if (map_parameters.skip_self || map_parameters.skip_prefix) {
            // before we set skch::filter::ONETOONE here
            // but this does not provide a clear benefit in all-to-all
            // as it sometimes introduces cases of over-filtering
            map_parameters.filterMode = skch::filter::MAP;
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

    if (!args::get(wfa_score_params).empty()) {
        const std::vector<std::string> params_str = skch::CommonFunc::split(args::get(wfa_score_params), ',');
        if (params_str.size() != 3) {
            std::cerr << "[wfmash] ERROR error: 3 scoring parameters must be given to -g/--wfa-params"//either 3 or 5 scoring parameters must be given to -g/--wflamda-params
                      << std::endl;
            exit(1);
        }

        std::vector<int> params(params_str.size());
        std::transform(params_str.begin(), params_str.end(), params.begin(),
                       [](const std::string &s) { return std::stoi(s); });

        align_parameters.wfa_mismatch_score = params[0];
        align_parameters.wfa_gap_opening_score = params[1];
        align_parameters.wfa_gap_extension_score = params[2];

        /*if (params.size() == 6) {
            align_parameters.wflambda_mismatch_score = params[0];
            align_parameters.wflambda_gap_opening_score = params[1];
            align_parameters.wflambda_gap_extension_score = params[2];
            xx = params[4];
            xx = params[5];
        }*/
    } else {
        align_parameters.wfa_mismatch_score = -1;
        align_parameters.wfa_gap_opening_score = -1;
        align_parameters.wfa_gap_extension_score = -1;
    }

    if (!args::get(wflign_score_params).empty()) {
        const std::vector<std::string> params_str = skch::CommonFunc::split(args::get(wflign_score_params), ',');
        if (params_str.size() != 3) {
            std::cerr << "[wfmash] ERROR error: 3 scoring parameters must be given to -G/--wflign-params."//either 3 or 5 scoring parameters must be given to -G/--wflign-params
                      << std::endl;
            exit(1);
        }

        std::vector<int> params(params_str.size());
        std::transform(params_str.begin(), params_str.end(), params.begin(),
                       [](const std::string &s) { return std::stoi(s); });

        align_parameters.wflign_mismatch_score = params[0];
        align_parameters.wflign_gap_opening_score = params[1];
        align_parameters.wflign_gap_extension_score = params[2];

        /*if (params.size() == 6) {
            align_parameters.wflign_gap_opening_score = params[0];
            align_parameters.wflign_gap_extension_score = params[1];
            align_parameters.wflign_gap_extension_score = params[2];
            xx = params[4];
            xx = params[5];
        }*/
    } else {
        align_parameters.wflign_mismatch_score = -1;
        align_parameters.wflign_gap_opening_score = -1;
        align_parameters.wflign_gap_extension_score = -1;
    }

    if (wflign_max_mash_dist) {
        if (args::get(wflign_max_mash_dist) <= 0 || args::get(wflign_max_mash_dist) > 1) {
            std::cerr << "[wfmash] ERROR, skch::parseandSave, max mash distance must be greater than 0 and less than or equal to 1." << std::endl;
            exit(1);
        }
        align_parameters.wflign_max_mash_dist = args::get(wflign_max_mash_dist);
    } else {
        align_parameters.wflign_max_mash_dist = -1;
    }

    if (wflign_min_wavefront_length) {
        if (args::get(wflign_min_wavefront_length) <= 0) {
            std::cerr << "[wfmash] ERROR, skch::parseandSave, min wavefront length for heuristic WFlign must be greater than 0." << std::endl;
            exit(1);
        }
        align_parameters.wflign_min_wavefront_length = args::get(wflign_min_wavefront_length);
    } else {
        align_parameters.wflign_min_wavefront_length = 1024;
    }

    if (wflign_max_distance_threshold) {
        if (args::get(wflign_max_distance_threshold) <= 0) {
            std::cerr << "[wfmash] ERROR, skch::parseandSave, max distance threshold for heuristic WFlign must be greater than 0." << std::endl;
            exit(1);
        }
        align_parameters.wflign_max_distance_threshold = args::get(wflign_max_distance_threshold);
    } else {
        align_parameters.wflign_max_distance_threshold = -1;
    }

    align_parameters.emit_md_tag = args::get(emit_md_tag);
    align_parameters.sam_format = args::get(sam_format);
    align_parameters.no_seq_in_sam = args::get(no_seq_in_sam);
    map_parameters.split = !args::get(no_split);
    map_parameters.dropRand = !args::get(keep_same_score);
    align_parameters.split = !args::get(no_split);

    map_parameters.mergeMappings = !args::get(no_merge);

    if (segment_length) {
        const int64_t s = wfmash::handy_parameter(args::get(segment_length));

        if (s <= 0) {
            std::cerr << "[wfmash] ERROR, skch::parseandSave, segment length has to be a float value greater than 0." << std::endl;
            exit(1);
        }

        if (s < 100) {
            std::cerr << "[wfmash] ERROR, skch::parseandSave, minimum segment length is required to be >= 100 bp." << std::endl
                      << "[wfmash] This is because Mashmap is not designed for computing short local alignments." << std::endl;
            exit(1);
        }
        map_parameters.segLength = s;
    } else {
        map_parameters.segLength = 5000;
    }

    if (map_pct_identity) {
        if (args::get(map_pct_identity) < 50) {
            std::cerr << "[wfmash] ERROR, skch::parseandSave, minimum nucleotide identity requirement should be >= 50\%." << std::endl;
            exit(1);
        }
        map_parameters.percentageIdentity = (float) (args::get(map_pct_identity)/100.0); // scale to [0,1]
    } else {
        map_parameters.percentageIdentity = skch::fixed::percentage_identity;
    }

    if (block_length) {
        const int64_t l = wfmash::handy_parameter(args::get(block_length));

        if (l < 0) {
            std::cerr << "[wfmash] ERROR, skch::parseandSave, min block length has to be a float value greater than or equal to 0." << std::endl;
            exit(1);
        }

        map_parameters.block_length = l;
    } else {
        // n.b. we map-merge across gaps up to 3x segment length
        // and then filter for things that are at least block_length long
        map_parameters.block_length = 5 * map_parameters.segLength;
    }

    if (chain_gap) {
        const int64_t l = wfmash::handy_parameter(args::get(chain_gap));
        if (l < 0) {
            std::cerr << "[wfmash] ERROR, skch::parseandSave, chain gap has to be a float value greater than or equal to 0." << std::endl;
            exit(1);
        }
        map_parameters.chain_gap = l;
        align_parameters.chain_gap = l;
    } else {
        map_parameters.chain_gap = 20000;
        align_parameters.chain_gap = 20000;
    }

    if (drop_low_map_pct_identity) {
        map_parameters.keep_low_pct_id = false;
    } else {
        map_parameters.keep_low_pct_id = true;
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
        map_parameters.kmerSize = 19;
    }

    if (kmer_pct_threshold) {
        map_parameters.kmer_pct_threshold = args::get(kmer_pct_threshold);
    } else {
        map_parameters.kmer_pct_threshold = 0.001; // in percent! so we keep 99.999% of kmers
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


    if (align_input_paf) {
        if (invert_filtering) {
            align_parameters.min_identity = map_parameters.percentageIdentity * 0.8; // in [0,1]
        } else {
            align_parameters.min_identity = 0; // disabled
        }
    } else {
        if (invert_filtering) {
            align_parameters.min_identity = 0; // disable
        } else {
            align_parameters.min_identity = map_parameters.percentageIdentity * 0.8; // in [0,1]
        }
    }

    if (wflambda_segment_length) {
        align_parameters.wflambda_segment_length = args::get(wflambda_segment_length);
    } else {
        align_parameters.wflambda_segment_length = 256;
    }

    if (wflign_max_len_major) {
        const uint64_t wflign_max_len_major_ = (uint64_t)wfmash::handy_parameter(args::get(wflign_max_len_major));

        if (wflign_max_len_major_ <= 0) {
            std::cerr << "[wfmash] ERROR, skch::parseandSave, maximum length to patch in the major axis has to be a float value greater than 0." << std::endl;
            exit(1);
        }

        align_parameters.wflign_max_len_major = wflign_max_len_major_;
    } else {
        align_parameters.wflign_max_len_major = map_parameters.segLength * 512;
    }

    if (wflign_max_len_minor) {
        const uint64_t wflign_max_len_minor_ = (uint64_t)wfmash::handy_parameter(args::get(wflign_max_len_minor));

        if (wflign_max_len_minor_ <= 0) {
            std::cerr << "[wfmash] ERROR, skch::parseandSave, maximum length to patch in the minor axis has to be a float value greater than 0." << std::endl;
            exit(1);
        }

        align_parameters.wflign_max_len_minor = wflign_max_len_minor_;
    } else {
        align_parameters.wflign_max_len_minor = map_parameters.segLength * 128;
    }

    if (wflign_erode_k) {
        align_parameters.wflign_erode_k = args::get(wflign_erode_k);
    } else {
        align_parameters.wflign_erode_k = -1; // will trigger estimation based on sequence divergence
    }

    if (wflign_max_patching_score) {
        align_parameters.wflign_max_patching_score = args::get(wflign_max_patching_score);
    } else {
        align_parameters.wflign_max_patching_score = 0; // will trigger estimation based on gap penalties and sequence length
    }

    if (thread_count) {
        map_parameters.threads = args::get(thread_count);
        align_parameters.threads = args::get(thread_count);
    } else {
        map_parameters.threads = 1;
        align_parameters.threads = 1;
    }

    // Compute optimal window size for sketching
    {
        const int64_t ss = sketch_size && args::get(sketch_size) >= 0 ? args::get(sketch_size) : -1;
        if (ss > 0) {
            map_parameters.sketchSize = ss;
        } else {
            const double md = 1 - map_parameters.percentageIdentity;
            double dens = 0.02 * (1 + (md / 0.05));
            map_parameters.sketchSize = dens * (map_parameters.segLength - map_parameters.kmerSize);
        }
    }

    if (kmer_complexity)
    {
        map_parameters.kmerComplexityThreshold = args::get(kmer_complexity);
    } else {
        map_parameters.kmerComplexityThreshold = 0;
    }

    map_parameters.filterLengthMismatches = true;

    map_parameters.stage1_topANI_filter = !bool(no_hg_filter); 

    if (hg_filter_ani_diff)
    {
        map_parameters.ANIDiff = args::get(hg_filter_ani_diff);
        map_parameters.ANIDiff /= 100;
    } else {
        map_parameters.ANIDiff = skch::fixed::ANIDiff;
    }

    if (hg_filter_conf)
    {
        map_parameters.ANIDiffConf = args::get(hg_filter_conf);
        map_parameters.ANIDiffConf /= 100;
    } else {
        map_parameters.ANIDiffConf = skch::fixed::ANIDiffConf;
    }

    //if (window_minimizers) {
        //map_parameters.world_minimizers = false;
    //} else {
        //map_parameters.world_minimizers = true;
    //}

    if (approx_mapping) {
        map_parameters.outFileName = "/dev/stdout";
        yeet_parameters.approx_mapping = true;
    } else {
        yeet_parameters.approx_mapping = false;

        if (tmp_base) {
            temp_file::set_dir(args::get(tmp_base));
        } else {
            char* cwd = get_current_dir_name();
            temp_file::set_dir(std::string(cwd));
            free(cwd);
        }

        if (align_input_paf) {
            yeet_parameters.remapping = true;
            map_parameters.outFileName = args::get(align_input_paf);
            align_parameters.mashmapPafFile = temp_file::create();
        } else {
            map_parameters.outFileName = temp_file::create();
            align_parameters.mashmapPafFile = map_parameters.outFileName;
        }
        align_parameters.pafOutputFile = "/dev/stdout";
    }

#ifdef WFA_PNG_AND_TSV
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

	map_parameters.legacy_output = false;

    //Check if files are valid
    skch::validateInputFiles(map_parameters.querySequences, map_parameters.refSequences);

    temp_file::set_keep_temp(args::get(keep_temp_files));

}

}
