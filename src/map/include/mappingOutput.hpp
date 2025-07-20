/**
 * @file    mappingOutput.hpp
 * @brief   Result processing and output for mappings
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef MAPPING_OUTPUT_HPP
#define MAPPING_OUTPUT_HPP

#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cmath>
#include <functional>
#include "map/include/base_types.hpp"
#include "map/include/sequenceIds.hpp"
#include "map/include/map_parameters.hpp"

namespace skch
{
  /**
   * @brief Class for handling mapping output and reporting
   */
  class MappingOutput {
  public:
    /**
     * @brief Check mapping boundaries and fix if needed
     */
    template <typename VecIn>
    static void mappingBoundarySanityCheck(InputSeqProgContainer* input, 
                                          VecIn &readMappings,
                                          const SequenceIdManager& idManager)
    {
      for(auto &e : readMappings)
      {
        //reference start pos
        {
          if(e.refStartPos < 0)
            e.refStartPos = 0;
          if(e.refStartPos >= idManager.getSequenceLength(e.refSeqId))
            e.refStartPos = idManager.getSequenceLength(e.refSeqId) - 1;
        }

        //reference end pos
        {
          if(e.refEndPos < e.refStartPos)
            e.refEndPos = e.refStartPos;
          if(e.refEndPos >= idManager.getSequenceLength(e.refSeqId))
            e.refEndPos = idManager.getSequenceLength(e.refSeqId) - 1;
        }

        //query start pos
        {
          if(e.queryStartPos < 0)
            e.queryStartPos = 0;
          if(e.queryStartPos >= input->len)
            e.queryStartPos = input->len;
        }

        //query end pos
        {
          if(e.queryEndPos < e.queryStartPos)
            e.queryEndPos = e.queryStartPos;
          if(e.queryEndPos >= input->len)
            e.queryEndPos = input->len;
        }
      }
    }

    /**
     * @brief Report final read mappings to output stream
     */
    static void reportReadMappings(MappingResultsVector_t &readMappings, 
                                  const std::string &queryName,
                                  std::ostream &outstrm,
                                  const SequenceIdManager& idManager,
                                  const Parameters& param,
                                  std::function<void(const MappingResult&)> processMappingResults = nullptr)
    {
      // Sort mappings by chain ID and query position
      std::sort(readMappings.begin(), readMappings.end(),
          [](const MappingResult &a, const MappingResult &b) {
              return std::tie(a.splitMappingId, a.queryStartPos)
                  < std::tie(b.splitMappingId, b.queryStartPos);
          });

      // Assign chain positions within each chain
      int current_chain = -1;
      int chain_pos = 0;
      int chain_length = 0;
      
      // First pass - count chain lengths
      for (size_t i = 0; i < readMappings.size(); ++i) {
          if (readMappings[i].splitMappingId != current_chain) {
              current_chain = readMappings[i].splitMappingId;
              chain_length = 1;
              // Count forward to find chain length
              for (size_t j = i + 1; j < readMappings.size(); ++j) {
                  if (readMappings[j].splitMappingId == current_chain) {
                      chain_length++;
                  } else {
                      break;
                  }
              }
              // Assign length to all mappings in this chain
              readMappings[i].chain_length = chain_length;
              chain_pos = 1;
          } else {
              readMappings[i].chain_length = chain_length;
              chain_pos++;
          }
          readMappings[i].chain_pos = chain_pos;
      }

      //Print the results
      for(auto &e : readMappings)
      {
        float fakeMapQ = e.nucIdentity == 1 ? 255 : std::round(-10.0 * std::log10(1-(e.nucIdentity)));
        std::string sep = param.legacy_output ? " " : "\t";

        outstrm  << (param.filterMode == filter::ONETOONE ? idManager.getSequenceName(e.querySeqId) : queryName)
                 << sep << e.queryLen
                 << sep << e.queryStartPos
                 << sep << e.queryEndPos - (param.legacy_output ? 1 : 0)
                 << sep << (e.strand == strnd::FWD ? "+" : "-")
                 << sep << idManager.getSequenceName(e.refSeqId)
                 << sep << idManager.getSequenceLength(e.refSeqId)
                 << sep << e.refStartPos
                 << sep << e.refEndPos - (param.legacy_output ? 1 : 0);

        if (!param.legacy_output) 
        {
          outstrm  << sep << e.conservedSketches
                   << sep << e.blockLength
                   << sep << fakeMapQ
                   << sep << "id:f:" << e.nucIdentity
                   << sep << "kc:f:" << e.kmerComplexity;
          if (!param.mergeMappings) 
          {
            outstrm << sep << "jc:f:" << float(e.conservedSketches) / e.sketchSize;
          } else {
            // Use chain_id if set, otherwise use splitMappingId
            int32_t chainId = (e.chain_id >= 0) ? e.chain_id : e.splitMappingId;
            outstrm << sep << "ch:Z:" << chainId << "." << e.chain_pos << "." << e.chain_length;
          }
        } else
        {
          outstrm << sep << e.nucIdentity * 100.0;
        }

#ifdef DEBUG
        outstrm << std::endl;
#else
        outstrm << "\n";
#endif

        //User defined processing of the results
        if(processMappingResults != nullptr)
          processMappingResults(e);
      }
    }

    /**
     * @brief Utility function to save L2 results to vector
     */
    static void insertL2ResultsToVec(MappingResultsVector_t &v, const MappingResult &reportedL2Result)
    {
      v.push_back(reportedL2Result);
    }
  };
}

#endif // MAPPING_OUTPUT_HPP