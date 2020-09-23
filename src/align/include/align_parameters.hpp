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
    float percentageIdentity;                     //user defined threshold for good similarity
    float min_identity;                           // drop alignments below this identity threshold
    int wf_min;                                   // minimum wavefront length to trigger WF_reduce wavefront pruning
    int wf_dist;                                  // max distance threshold that a wavefront may lag behind the best wavefront and not be removed

    std::vector<std::string> refSequences;        //reference sequence(s)
    std::vector<std::string> querySequences;      //query sequence(s)
    std::string mashmapPafFile;                   //mashmap paf mapping file
    std::string pafOutputFile;                    //sam output file name
};

}

#endif
