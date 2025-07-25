/**
 * @file    map_parameters.hpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef SKETCH_CONFIG_HPP
#define SKETCH_CONFIG_HPP

#include <vector>
#include <unordered_set>
#include <filesystem>
namespace stdfs = std::filesystem;

#include "common/ALeS.hpp"
#include "base_types.hpp"

namespace skch
{


struct ales_params {
  uint32_t weight{} ;
  uint32_t seed_count{};
  float similarity{};
  uint32_t region_length{};
};

/**
 * @brief   configuration parameters for building sketch
 *          expected to be initialized using command line arguments
 */
struct Parameters
{
    int kmerSize;                                     //kmer size for sketching
    offset_t windowLength;                             //window size for sketching
    offset_t block_length;                             // minimum (potentially merged) block to keep if we aren't split
    offset_t chain_gap;                                // max distance for 2d range union-find mapping chaining
    uint64_t max_mapping_length;                      // maximum length of a mapping
    int alphabetSize;                                 //alphabet size
    offset_t referenceSize;                           //Approximate reference size
    float percentageIdentity;                         //user defined threshold for good similarity
    bool stage2_full_scan;                            //Instead of using the best intersection for a given candidate region, compute the minhash for every position in the window
    bool stage1_topANI_filter;                        //Use the ANI filter in stage 1
    float ANIDiff;                                    //ANI distance threshold below best mapping to retain in stage 1 filtering
    float ANIDiffConf;                                //Confidence of stage 1 ANI filtering threshold
    int filterMode;                                   //filtering mode in mashmap
    uint32_t numMappingsForSegment;                   //how many mappings to retain for each segment
    uint32_t numMappingsForShortSequence;             //how many secondary alignments we keep for reads < segLength
    bool dropRand;                                    //drop mappings w/ same score until only numMappingsForSegment remain
    int threads;                                      //execution thread count
    std::vector<std::string> refSequences;            //reference sequence(s)
    std::vector<std::string> querySequences;          //query sequence(s)
    std::string outFileName;                          //output file name
    stdfs::path indexFilename;                        //output file name of index
    bool overwrite_index;                             //overwrite index if it exists
    bool create_index_only;                           //only create index and exit
    bool split;                                       //Split read mapping (done if this is true)
    bool lower_triangular;                            // set to true if we should filter out half of the mappings
    bool skip_self;                                   //skip self mappings
    bool skip_prefix;                                 //skip mappings to sequences with the same prefix
    char prefix_delim;                                //the prefix delimiter
    std::string target_list;                          //file containing list of target sequences
    std::string target_prefix;                        //prefix for target sequences to use
    bool mergeMappings;                               //if we should merge consecutive segment mappings
    bool keep_low_pct_id;                             //true if we should keep mappings whose estimated identity < percentageIdentity
    bool report_ANI_percentage;                       //true if ANI should be in [0,100] as opposed to [0,1] (this is necessary for wfmash
    bool filterLengthMismatches;                      //true if filtering out length mismatches
    float kmerComplexityThreshold;                    //minimum kmer complexity to consider (default 0)

    std::string query_list;                           // file containing list of query sequence names
    std::vector<std::string> query_prefix;            // prefix for query sequences to use

    int sketchSize;
    double hgNumerator = 1.0;                         // Numerator for the hypergeometric filter's Jaccard similarity calculation
    uint64_t totalReferenceSize = 0;                  // Total size of all reference sequences
    uint64_t estimatedUniqueKmers = 0;                // Estimate of total unique k-mers
    bool use_spaced_seeds;                            //
    ales_params spaced_seed_params;                   //
    double spaced_seed_sensitivity;                   //
    std::vector<ales::spaced_seed> spaced_seeds;      //
    bool world_minimizers;
    uint64_t sparsity_hash_threshold;                 // keep mappings that hash to <= this value
    double overlap_threshold;                         // minimum overlap for a mapping to be considered
    int64_t scaffold_max_deviation;                  // max diagonal deviation from scaffold chains
    int64_t scaffold_gap;                           // gap threshold for scaffold chaining
    int64_t scaffold_min_length = 50000;            // minimum scaffold block length
    std::string scaffold_output_file;               // optional file to output scaffold mappings
    
    bool legacy_output;
    //std::unordered_set<std::string> high_freq_kmers;  //
    int64_t index_by_size = std::numeric_limits<int64_t>::max();  // Target total size of sequences for each index subset
    int minimum_hits = -1;  // Minimum number of hits required for L1 filtering (-1 means auto)
    double max_kmer_freq = 0.0002;  // Maximum allowed k-mer frequency fraction (0-1) or count (>1)

    bool use_progress_bar = false;
    bool auto_pct_identity = true;  // default to auto identity estimation
    int ani_percentile = 25;  // which percentile to use (25, 50, 75, etc.)
    float ani_adjustment = -5.0;  // adjustment to apply to the percentile (+/- percentage points)
    bool use_streaming_minhash = true;  // use efficient streaming MinHash algorithm (default enabled)
    int ani_sketch_size = 1000;  // sketch size for ANI estimation
};


/**
 * @brief     Default values or internal figures not exposed at the command line interface
 */
namespace fixed
{

//float filter_score_best_range = .99;              // mapping score above a certain fraction of best score is
//considered good by filtering algorithm

//int max_best_mappings_per_position = 25;          // At a particular position, if algorithm finds more than a certain best
//mappings, it doesn't mark them as best anymore

double ss_table_max = 1000.0;                       // Maximum size of dp table for filtering
double pval_cutoff = 1e-3;                          // p-value cutoff for determining window size
float confidence_interval = 0.95;                   // Confidence interval to relax jaccard cutoff for mapping (0-1)
float percentage_identity = 0.70;                   // Percent identity in the mapping step
float ANIDiff = 0.0;                                // Stage 1 ANI diff threshold
float ANIDiffConf = 0.999;                          // ANI diff confidence
std::string VERSION = "3.5.0";                      // Version of MashMap
}
}

#endif
