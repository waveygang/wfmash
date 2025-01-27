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
    uint16_t rankMapping;             //rank of the mapping for the query qId. It has the same variable type of num_mappings_for_segments (used for the SAM output format)

    std::string qId;                    //query sequence(s) 
    std::string refId;                  //reference sequence(s)
    skch::offset_t qStartPos;           //mapping boundary start offset on query
    skch::offset_t qEndPos;             //mapping boundary end offset on query
    skch::offset_t rStartPos;           //mapping boundary start offset on ref
    skch::offset_t rEndPos;             //mapping boundary end offset on ref
    skch::strand_t strand;              //mapping strand
    float mashmap_estimated_identity;

    // Chain metadata
    int32_t chain_id{-1};               // Unique ID for this chain (-1 if not part of chain)
    int32_t chain_length{1};            // Total segments in chain (1 if not part of chain)
    int32_t chain_pos{1};               // Position in chain, 1-based (1 if not part of chain)
  };

  typedef std::unordered_map <std::string, std::string> refSequenceMap_t;
}

#endif
