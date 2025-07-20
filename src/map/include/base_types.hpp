/**
 * @file    base_types.hpp
 * @brief   Critical type definitions for mapping algorithm
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef BASE_TYPES_MAP_HPP
#define BASE_TYPES_MAP_HPP

#include <tuple>
#include <vector>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <limits>
#include "common/progress.hpp"

namespace skch
{
  // Forward declaration
  class SequenceIdManager;
  typedef uint64_t hash_t;    //hash type
  typedef int64_t offset_t;   //position within sequence
  typedef int32_t seqno_t;    //sequence counter in file
  typedef int16_t strand_t;   //sequence strand 
  typedef int8_t side_t;      //sequence strand 

  //C++ timer
  typedef std::chrono::high_resolution_clock Time;

  //Information about each minimizer
  struct MinmerInfo
  {
    hash_t hash;                              //hash value
    offset_t wpos;                            //First (left-most) window position when the minimizer is saved
    offset_t wpos_end;
    seqno_t seqId;                            //sequence or contig id
    strand_t strand;                          //strand information

    //Lexographical equality comparison
    bool operator ==(const MinmerInfo& x) const {
      return std::tie(hash, seqId, wpos, strand) 
        == std::tie(x.hash, x.seqId, x.wpos, x.strand);
    }

    bool operator !=(const MinmerInfo& x) const {
      return std::tie(hash, seqId, wpos, strand) 
        != std::tie(x.hash, x.seqId, x.wpos, x.strand);
    }

    static bool equalityByHash(const MinmerInfo& x, const MinmerInfo& y) {
      return x.hash == y.hash;
    }

    static bool lessByHash(const MinmerInfo& x, const MinmerInfo& y) {
      return x.hash < y.hash;
    }

    // Sort based on start point
    bool operator <(const MinmerInfo& x) const {
      return std::tie(seqId, wpos) 
        < std::tie(x.seqId, x.wpos);
    }
  };

  // Endpoints for minmer intervals
  struct IntervalPoint
  {
    offset_t pos;
    hash_t hash;
    seqno_t seqId;
    side_t side;

    // Sort interval points. 
    // For a pair of points at the same seqId/pos, the end point should be first
    bool operator <(const IntervalPoint& x) const {
      return std::tie(seqId, pos, side) 
        < std::tie(x.seqId, x.pos, x.side);
    }
  };

  template <class It>
  struct boundPtr {
    It it;
    It end;

    bool operator<(const boundPtr& other) const {
      return *it < *(other.it);
    }
  };


  typedef hash_t MinmerMapKeyType;
  typedef std::vector<IntervalPoint> MinmerMapValueType;

  //Metadata recording for contigs in the reference DB
  struct ContigInfo
  {
    std::string name;       //Name of the sequence
    offset_t len;           //Length of the sequence
    int groupId;            //Group ID for the sequence
  };

  //Label tags for strand information
  enum strnd : strand_t
  {
    FWD = 1,
    AMBIG = 0,
    REV = -1
  };

  enum event : int
  {
    BEGIN = 1,
    END = 2
  };

  //filter mode in mashmap
  enum filter : int
  {
    MAP = 1,                              //filter by query axis
    ONETOONE = 2,                         //filter by query axis and reference axis
    NONE = 3                              //no filtering
  };

  // Enum for tracking which side of an interval a point represents
  enum side : side_t
  {
    OPEN = 1,  
    CLOSE = -1
  };  

  struct SeqCoord
  {
    seqno_t seqId;
    offset_t pos;
  };

  struct KmerInfo
  {
    hash_t hash;
    seqno_t seqId;
    offset_t pos;
    strand_t strand; 
  };

  template <class T>
  inline void hash_combine(std::size_t & s, const T & v)
  {
      std::hash<T> h;
      s^= h(v) + 0x9e3779b9 + (s<< 6) + (s>> 2);
  }

  //Fragment mapping result
  //Do not save variable sized objects in this struct
  struct MappingResult
  {
    offset_t queryLen;                                  //length of the query sequence
    offset_t refStartPos;                               //start position of the mapping on reference
    offset_t refEndPos;                                 //end pos
    offset_t queryStartPos;                             //start position of the query for this mapping
    offset_t queryEndPos;                               //end position of the query for this mapping
    seqno_t refSeqId;                                   //internal sequence id of the reference contig
    seqno_t querySeqId;                                 //internal sequence id of the query sequence
    offset_t blockLength;                                    //the block length of the mapping
    float blockNucIdentity;
          
    float nucIdentity;                                  //calculated identity
    float nucIdentityUpperBound;                        //upper bound on identity (90% C.I.)
    int sketchSize;                                     //sketch size
    int conservedSketches;                              //count of conserved sketches
    strand_t strand;                                    //strand
    int approxMatches;                                  //the approximate number of matches in the alignment

                                                        //--for split read mapping

    long double kmerComplexity;                               // Estimated sequence complexity
    int n_merged;                                       // how many mappings we've merged into this one
    offset_t splitMappingId;                            // To identify split mappings that are chained
    uint8_t discard;                                    // set to 1 for deletion
    bool overlapped;                                    // set to true if this mapping is overlapped with another mapping
    bool selfMapFilter;                                 // set to true if a long-to-short mapping in all-vs-all mode (we report short as the query)
    double chainPairScore;                              // best score for potential chain pair
    int64_t chainPairId;                                // best partner mapping for potential chain pair
    int32_t chain_id{-1};                               //unique ID for this chain (-1 if not part of chain)
    int32_t chain_length{1};                            //total segments in chain (1 if not part of chain)
    int32_t chain_pos{1};                               //position in chain, 1-based (1 if not part of chain)

    offset_t qlen() {                                   //length of this mapping on query axis
      return queryEndPos - queryStartPos + 1;
    }

    offset_t rlen() {                                   //length of this mapping on reference axis
      return refEndPos - refStartPos + 1;
    }

    size_t hash(void) {
      size_t res = 0;
      hash_combine(res, queryLen);
      hash_combine(res, refStartPos);
      hash_combine(res, refEndPos);
      hash_combine(res, queryStartPos);
      hash_combine(res, queryEndPos);
      hash_combine(res, refSeqId);
      hash_combine(res, querySeqId);
      hash_combine(res, blockLength);
      hash_combine(res, nucIdentity);
      hash_combine(res, nucIdentityUpperBound);
      hash_combine(res, sketchSize);
      hash_combine(res, conservedSketches);
      hash_combine(res, strand);
      hash_combine(res, approxMatches);
      return res;
    }

  };

  typedef std::vector<MappingResult> MappingResultsVector_t;

  // Minimal mapping representation for memory efficiency
  // Stores only essential data in 32 bytes instead of 176+ bytes
  struct MinimalMapping {
    uint32_t ref_seqId;     // 4 bytes - reference sequence ID
    uint32_t ref_pos;       // 4 bytes - reference start position
    uint32_t query_pos;     // 4 bytes - query start position  
    uint32_t n_merged;      // 4 bytes - number of merged mappings (was uint16_t, now uint32_t for large scaffolds)
    uint32_t conservedSketches; // 4 bytes - conserved sketches count (was uint16_t, now uint32_t for large sequences)
    uint32_t length;        // 4 bytes - mapping length (was uint16_t, now uint32_t for long reads)
    uint16_t identity;      // 2 bytes - identity percentage (0-10000 for 0.00-100.00%)
    uint8_t flags;          // 1 byte - packed flags (strand, discard, overlapped)
    uint8_t kmerComplexity; // 1 byte - kmer complexity (0-255)
                           // Total: 32 bytes exactly, no padding needed
    
    // Flag bit positions
    static constexpr uint8_t FLAG_STRAND_MASK = 0x03;  // bits 0-1 for strand (-1, 0, 1)
    static constexpr uint8_t FLAG_DISCARD = 0x04;      // bit 2
    static constexpr uint8_t FLAG_OVERLAPPED = 0x08;   // bit 3
    
    // Helper methods for flag manipulation
    strand_t getStrand() const {
      uint8_t strand_bits = flags & FLAG_STRAND_MASK;
      return (strand_bits == 0) ? strnd::REV : 
             (strand_bits == 1) ? strnd::AMBIG : strnd::FWD;
    }
    
    void setStrand(strand_t s) {
      flags = (flags & ~FLAG_STRAND_MASK) | 
              ((s == strnd::REV) ? 0 : (s == strnd::AMBIG) ? 1 : 2);
    }
    
    bool isDiscarded() const { return flags & FLAG_DISCARD; }
    void setDiscarded(bool d) { 
      if (d) flags |= FLAG_DISCARD; 
      else flags &= ~FLAG_DISCARD;
    }
    
    bool isOverlapped() const { return flags & FLAG_OVERLAPPED; }
    void setOverlapped(bool o) { 
      if (o) flags |= FLAG_OVERLAPPED;
      else flags &= ~FLAG_OVERLAPPED;
    }
  };

  // Forward declarations for conversion functions
  MappingResult expandMinimalMapping(const MinimalMapping& m, 
                                   seqno_t querySeqId, 
                                   offset_t queryLen);
  
  MinimalMapping compressMapping(const MappingResult& full);

  //Vector type for storing MinmerInfo
  typedef std::vector<MinmerInfo> MinVec_Type;

  //Container to save copy of kseq object
  struct InputSeqContainer
  {
    seqno_t seqId;                              //sequence id
    offset_t len;                               //sequence length
    std::string seq;                            //sequence string
    std::string name;                        //sequence name


    /*
     * @brief               constructor
     * @param[in] kseq_seq  complete read or reference sequence
     * @param[in] kseq_id   sequence id name
     * @param[in] len       length of sequence
     */
      InputSeqContainer(const std::string& s, const std::string& name, seqno_t id)
          : seqId(id)
          , len(s.length())
          , seq(s)
          , name(name) { }
  };

  struct InputSeqProgContainer : InputSeqContainer
  {
    using InputSeqContainer::InputSeqContainer;
    progress_meter::ProgressMeter& progress;    //progress meter (shared)
                                                

    /*
     * @brief               constructor
     * @param[in] kseq_seq  complete read or reference sequence
     * @param[in] kseq_id   sequence id name
     * @param[in] len       length of sequence
     */
      InputSeqProgContainer(const std::string& s, const std::string& name, seqno_t id, progress_meter::ProgressMeter& pm)
          : InputSeqContainer(s, name, id)
          , progress(pm) { }
  };


  //Output type of map function
  struct MapModuleOutput
  {
    MappingResultsVector_t readMappings;  //read mapping coordinates
    std::string qseqName;                 //query sequence id
    offset_t qseqLen;                     //query sequence length

    //Function to erase all output mappings
    void reset()
    {
      this->readMappings.clear();
    }
  };

  //Information about fragment sequence during L1/L2 mapping
  template <typename MinmerVec>
    struct QueryMetaData
    {
      char *seq;                          //query sequence pointer
      seqno_t seqId;                      //query sequence id
      offset_t len;                       //length of this query sequence
      offset_t fullLen;                   //length of the full sequence it derives from
      int sketchSize;                     //sketch size
      std::string seqName;                //sequence name
      MinmerVec minmerTableQuery;         //Vector of minmers in the query
      MinmerVec seedHits;                 //Vector of minmers in the reference
      int refGroup;                       //Prefix group of sequence
      float kmerComplexity;                //Estimated sequence complexity
    };

  // Implementation of conversion functions
  inline MinimalMapping compressMapping(const MappingResult& full) {
    MinimalMapping minimal;
    
    // Direct assignments for simple fields
    minimal.ref_seqId = static_cast<uint32_t>(full.refSeqId);
    minimal.ref_pos = static_cast<uint32_t>(full.refStartPos);
    minimal.query_pos = static_cast<uint32_t>(full.queryStartPos);
    
    // Calculate length (use the larger of ref/query length)
    offset_t ref_len = full.refEndPos - full.refStartPos;
    offset_t query_len = full.queryEndPos - full.queryStartPos;
    minimal.length = static_cast<uint32_t>(std::max(ref_len, query_len));
    
    // Store n_merged (no cap needed with uint32_t)
    minimal.n_merged = full.n_merged;
    
    // Convert identity from float (0.0-1.0) to uint16_t (0-10000 for 0.00-100.00%)
    minimal.identity = static_cast<uint16_t>(std::round(full.nucIdentity * 10000));
    
    // Store conserved sketches (no cap needed with uint32_t)
    minimal.conservedSketches = full.conservedSketches;
    
    // Convert kmer complexity - can be > 1.0, so cap at 255
    minimal.kmerComplexity = static_cast<uint8_t>(std::min(255, static_cast<int>(std::round(full.kmerComplexity * 100))));
    
    // Pack flags
    minimal.flags = 0;
    minimal.setStrand(full.strand);
    minimal.setDiscarded(full.discard != 0);
    minimal.setOverlapped(full.overlapped);
    
    
    return minimal;
  }
  
  inline MappingResult expandMinimalMapping(const MinimalMapping& m, 
                                          seqno_t querySeqId, 
                                          offset_t queryLen) {
    MappingResult full;
    
    // Query information
    full.querySeqId = querySeqId;
    full.queryLen = queryLen;
    full.queryStartPos = m.query_pos;
    full.queryEndPos = m.query_pos + m.length;
    
    // Reference information
    full.refSeqId = m.ref_seqId;
    full.refStartPos = m.ref_pos;
    full.refEndPos = m.ref_pos + m.length;
    
    // Identity (convert from 0-10000 to 0.0-1.0)
    full.nucIdentity = m.identity / 10000.0f;
    full.nucIdentityUpperBound = full.nucIdentity; // Conservative estimate
    
    // Block information
    full.blockLength = m.length;
    full.blockNucIdentity = full.nucIdentity;
    full.approxMatches = static_cast<int>(std::round(full.nucIdentity * full.blockLength));
    
    // Strand and flags
    full.strand = m.getStrand();
    full.discard = m.isDiscarded() ? 1 : 0;
    full.overlapped = m.isOverlapped();
    
    // Restore metadata from MinimalMapping
    full.sketchSize = 0; // Not stored, would need to be recomputed
    full.conservedSketches = m.conservedSketches;
    full.kmerComplexity = m.kmerComplexity / 100.0f;
    full.n_merged = m.n_merged;
    full.splitMappingId = 0;
    full.selfMapFilter = false;
    full.chainPairScore = std::numeric_limits<double>::max();
    full.chainPairId = std::numeric_limits<int64_t>::min();
    full.chain_id = -1;
    full.chain_length = 1;
    full.chain_pos = 1;
    
    return full;
  }
}

#endif
