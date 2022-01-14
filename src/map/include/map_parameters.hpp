/**
 * @file    map_parameters.hpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef SKETCH_CONFIG_HPP 
#define SKETCH_CONFIG_HPP

#include <vector>
#include <unordered_set>

#include "common/ALeS.hpp"

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
    int64_t windowSize;                               //window size used for sketching
    int64_t segLength;                                //For split mapping case, this represents the fragment length
                                                      //for noSplit, it represents minimum read length to multimap
    int64_t block_length_min;                         // minimum (potentially merged) block to keep if we aren't split
    int alphabetSize;                                 //alphabet size
    uint64_t referenceSize;                           //Approximate reference size
    float percentageIdentity;                         //user defined threshold for good similarity
    int filterMode;                                   //filtering mode in mashmap
    uint16_t numMappingsForSegment;                   //how many mappings to retain for each segment
    uint16_t numMappingsForShortSequence;             //how many secondary alignments we keep for reads < segLength
    int threads;                                      //execution thread count
    std::vector<std::string> refSequences;            //reference sequence(s)
    std::vector<std::string> querySequences;          //query sequence(s)
    std::string outFileName;                          //output file name
    bool split;                                       //Split read mapping (done if this is true)
    bool skip_self;                                   //skip self mappings
    bool skip_prefix;                                 //skip mappings to sequences with the same prefix
    char prefix_delim;                                //the prefix delimiter
    bool mergeMappings;                               //if we should merge consecutive segment mappings
    bool keep_low_pct_id;                             //true if we should keep mappings whose estimated identity < percentageIdentity

    bool use_spaced_seeds;                            //
    ales_params spaced_seed_params;                   //
    double spaced_seed_sensitivity;                   //
    std::vector<ales::spaced_seed> spaced_seeds;      //

    //std::unordered_set<std::string> high_freq_kmers;  //
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

double pval_cutoff = 0.0;//1e-120;                  // p-value cutoff for determining window size
float confidence_interval = 0.95;                   // Confidence interval to relax jaccard cutoff for mapping (0-1)
float percentage_identity = 0.95;                   // Percent identity in the mapping step
}
}

#endif
