/**
 * @file    align_parameters.hpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef ALIGN_PARAMETERS_HPP
#define ALIGN_PARAMETERS_HPP

#include <vector>

namespace align {
  /**
   * @brief   parameters for generating mashmap alignments
   */

struct Parameters {
    int threads;                                  //execution thread count
    //float percentageIdentity;                     //user defined threshold for good similarity
    float min_identity;                           // drop alignments below this identity threshold
    uint64_t min_alignment_length;                // drop alignments shorter than this length (in bp)
    float min_block_identity;                     // drop alignments with block identity below this threshold
    //int wf_min;                                   // minimum wavefront length to trigger WF_reduce wavefront pruning
    //int wf_diff;                                  // max distance threshold that a wavefront may lag behind the best wavefront and not be removed
    //bool exact_wfa;                               // use exact WFA, avoiding adaptive wavefront reduction
    bool split;                                       //Split read mapping (done if this is true)

    //wflambda
    uint16_t wflambda_segment_length;             //segment length for wflambda

    bool force_biwfa_alignment;				   //force biwfa alignment
    bool force_wflign;                          //force alignment with WFlign instead of the default biWFA

    int wfa_mismatch_score;
    int wfa_gap_opening_score;
    int wfa_gap_extension_score;

    int wfa_patching_mismatch_score;
    int wfa_patching_gap_opening_score1;
    int wfa_patching_gap_extension_score1;
    int wfa_patching_gap_opening_score2;
    int wfa_patching_gap_extension_score2;

    // wflign
    int wflign_mismatch_score;
    int wflign_gap_opening_score;
    int wflign_gap_extension_score;
    float wflign_max_mash_dist;
    int wflign_min_wavefront_length;
    int wflign_max_distance_threshold;
    uint64_t wflign_max_len_major;
    uint64_t wflign_max_len_minor;
    int wflign_erode_k;
    int kmerSize;                                 //kmer size for pre-checking before aligning a fragment
    int64_t chain_gap;                            //max distance for 2d range union-find mapping chaining;
    int wflign_min_inv_patch_len;                 //minimum length of an inverted patch
    int wflign_max_patching_score;                //maximum score allowed for patching

    std::vector<std::string> refSequences;        //reference sequence(s)
    std::vector<std::string> querySequences;      //query sequence(s)
    std::string mashmapPafFile;                   //mashmap paf mapping file
    std::string pafOutputFile;                    //paf/sam output file name

    bool emit_md_tag;                             //Output the MD tag
    bool sam_format;                              //Emit the output in SAM format (PAF default)
    bool no_seq_in_sam;                           //Do not fill the SEQ field in SAM format
    bool disable_chain_patching;                  //Disable alignment patching at chain boundaries
    bool multithread_fasta_input;                 //Multithreaded fasta input
    uint64_t target_padding;                      //Additional padding around target sequence
    uint64_t query_padding;                       //Additional padding around query sequence

    bool use_progress_bar = false;

#ifdef WFA_PNG_TSV_TIMING
    // plotting
    std::string tsvOutputPrefix;                  //tsv files with wavefront information for each alignment
    uint64_t wfplot_max_size;                     // Max size (width/height) of the wfplot
    std::string prefix_wavefront_plot_in_png;     // Prefix of PNG files with wavefront plot for each alignment

    std::string path_patching_info_in_tsv;        // TSV file with patching statistics
#endif
};

}

#endif
