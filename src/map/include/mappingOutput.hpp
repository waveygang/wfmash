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

        //reference end pos - calculated from blockLength
        {
          if(e.refEndPos() < e.refStartPos)
            e.blockLength = 0; // Set blockLength to 0 if invalid
          if(e.refEndPos() >= idManager.getSequenceLength(e.refSeqId))
            e.blockLength = idManager.getSequenceLength(e.refSeqId) - 1 - e.refStartPos;
        }

        //query start pos
        {
          if(e.queryStartPos < 0)
            e.queryStartPos = 0;
          if(e.queryStartPos >= input->len)
            e.queryStartPos = input->len;
        }

        //query end pos - calculated from blockLength
        {
          if(e.queryEndPos() < e.queryStartPos)
            e.blockLength = 0; // Set blockLength to 0 if invalid
          if(e.queryEndPos() >= input->len)
            e.blockLength = input->len - e.queryStartPos;
        }
      }
    }

    /**
     * @brief Report final read mappings to output stream with chain information
     */
    static void reportReadMappings(MappingResultsVector_t &readMappings,
                                  const ChainInfoVector_t &chainInfo,
                                  const std::string &queryName,
                                  std::ostream &outstrm,
                                  const SequenceIdManager& idManager,
                                  const Parameters& param,
                                  std::function<void(const MappingResult&)> processMappingResults = nullptr,
                                  offset_t queryLen = 0)
    {
      // Sort mappings and chain info together by query position
      std::vector<size_t> indices(readMappings.size());
      std::iota(indices.begin(), indices.end(), 0);
      std::sort(indices.begin(), indices.end(),
          [&readMappings](size_t a, size_t b) {
              return readMappings[a].queryStartPos < readMappings[b].queryStartPos;
          });

      // Print the results
      for(size_t idx : indices)
      {
        auto &e = readMappings[idx];
        auto &chain = chainInfo[idx];
        
        float fakeMapQ = e.getNucIdentity() == 1 ? 255 : std::round(-10.0 * std::log10(1-(e.getNucIdentity())));
        std::string sep = param.legacy_output ? " " : "\t";

        outstrm  << queryName
                 << sep << queryLen
                 << sep << e.queryStartPos
                 << sep << e.queryEndPos() - (param.legacy_output ? 1 : 0)
                 << sep << (e.strand() == strnd::FWD ? "+" : "-")
                 << sep << idManager.getSequenceName(e.refSeqId)
                 << sep << idManager.getSequenceLength(e.refSeqId)
                 << sep << e.refStartPos
                 << sep << e.refEndPos() - (param.legacy_output ? 1 : 0);

        if (!param.legacy_output) 
        {
          outstrm  << sep << e.conservedSketches
                   << sep << e.blockLength
                   << sep << fakeMapQ
                   << sep << "id:f:" << e.getNucIdentity()
                   << sep << "kc:f:" << e.getKmerComplexity();
          if (!param.mergeMappings) 
          {
            outstrm << sep << "jc:f:" << 0.0;
          } else {
            outstrm << sep << "ch:Z:" << chain.chainId << "." << chain.chainPos << "." << chain.chainLen;
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

        // User defined processing of the results
        if(processMappingResults != nullptr)
          processMappingResults(e);
      }
    }

    /**
     * @brief Report final read mappings to output stream (legacy interface without chain info)
     */
    static void reportReadMappings(MappingResultsVector_t &readMappings, 
                                  const std::string &queryName,
                                  std::ostream &outstrm,
                                  const SequenceIdManager& idManager,
                                  const Parameters& param,
                                  std::function<void(const MappingResult&)> processMappingResults = nullptr,
                                  offset_t queryLen = 0)
    {
      // Create default chain info - each mapping is its own chain
      ChainInfoVector_t defaultChainInfo(readMappings.size());
      for (size_t i = 0; i < readMappings.size(); ++i) {
        defaultChainInfo[i] = {static_cast<uint32_t>(i), 1, 1};
      }
      
      // Call the version with chain info
      reportReadMappings(readMappings, defaultChainInfo, queryName, outstrm, 
                        idManager, param, processMappingResults, queryLen);
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