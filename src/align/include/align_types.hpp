/**
 * @file    align_types.hpp
 * @brief   Critical type defintions for generating alignments
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef ALIGN_TYPES_MAP_HPP 
#define ALIGN_TYPES_MAP_HPP

#include <tuple>
#include <unordered_map>

namespace align
{
  //Type for map value type used for
  //L1 stage lookup index
  struct MappingBoundaryRow
  {
    std::string qId;                    //query sequence(s)
    std::string refId;                  //reference sequence(s)
    skch::offset_t qStartPos;           //mapping boundary start offset on query
    skch::offset_t qEndPos;             //mapping boundary end offset on query
    skch::offset_t rStartPos;           //mapping boundary start offset on ref
    skch::offset_t rEndPos;             //mapping boundary end offset on ref
    skch::strand_t strand;              //mapping strand
    float mashmap_estimated_identity;

    // group_mapping and rank_mapping to identify uniquely each mapping
    uint64_t group_mapping;             //Group of mappings that refer to the same query portion
    uint32_t rank_mapping;              //Rank of the mapping in a group. It has the same variable type of num_mappings_for_segments (used for the SAM output format)
  };

  typedef std::unordered_map <std::string, std::string> refSequenceMap_t;
}

#endif
