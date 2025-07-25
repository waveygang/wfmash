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
#include "common/progress.hpp"

namespace skch
{
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

  /**
   * @struct MappingResult
   * @brief Compact 32-byte mapping structure for memory efficiency
   */
  struct MappingResult
  {
    uint32_t refSeqId;         // 4 bytes - reference sequence ID
    uint32_t refStartPos;      // 4 bytes - reference start position  
    uint32_t queryStartPos;    // 4 bytes - query start position
    uint32_t blockLength;      // 4 bytes - mapping block length
    uint32_t n_merged;         // 4 bytes - number of merged segments
    uint32_t conservedSketches;// 4 bytes - count of conserved sketches
    uint16_t nucIdentity;      // 2 bytes - scaled identity (0-10000 for 0.00-100.00%)
    uint8_t  flags;            // 1 byte - bit-packed flags (strand, discard, overlapped)
    uint8_t  kmerComplexity;   // 1 byte - scaled kmer complexity (0-100)
                               // Total: 32 bytes

    // Helper methods to extract flag values
    strand_t strand() const {
      return (flags & 0x01) ? strnd::REV : strnd::FWD;
    }
    
    bool discard() const {
      return (flags & 0x02) != 0;
    }
    
    bool overlapped() const {
      return (flags & 0x04) != 0;
    }
    
    // Helper methods to set flag values
    void setStrand(strand_t s) {
      if (s == strnd::REV) flags |= 0x01;
      else flags &= ~0x01;
    }
    
    void setDiscard(bool d) {
      if (d) flags |= 0x02;
      else flags &= ~0x02;
    }
    
    void setOverlapped(bool o) {
      if (o) flags |= 0x04;
      else flags &= ~0x04;
    }
    
    // Helper methods to get unscaled values
    float getNucIdentity() const {
      return nucIdentity / 10000.0f;
    }
    
    float getKmerComplexity() const {
      return kmerComplexity / 100.0f;
    }
    
    // Helper methods to set scaled values
    void setNucIdentity(float identity) {
      nucIdentity = static_cast<uint16_t>(roundf(identity * 10000.0f));
    }
    
    void setKmerComplexity(float complexity) {
      kmerComplexity = static_cast<uint8_t>(roundf(complexity * 100.0f));
    }

    // Calculate derived values
    offset_t refEndPos() const {
      return refStartPos + blockLength;
    }
    
    offset_t queryEndPos() const {
      return queryStartPos + blockLength;
    }
    
    offset_t qlen() const {
      return blockLength;
    }
    
    offset_t rlen() const {
      return blockLength;
    }
    
    // Compatibility accessors for code that expects these as fields
    float blockNucIdentity() const {
      return getNucIdentity();
    }
    
    // Hash function for compatibility
    size_t hash(void) const {
      size_t res = 0;
      hash_combine(res, refSeqId);
      hash_combine(res, refStartPos);
      hash_combine(res, queryStartPos);
      hash_combine(res, blockLength);
      hash_combine(res, nucIdentity);
      hash_combine(res, conservedSketches);
      hash_combine(res, flags);
      return res;
    }

    // Initialize all fields to default values
    MappingResult() : refSeqId(0), refStartPos(0), queryStartPos(0), 
                      blockLength(0), n_merged(1), conservedSketches(0),
                      nucIdentity(0), flags(0), kmerComplexity(0) {}
  };

  typedef std::vector<MappingResult> MappingResultsVector_t;

  /**
   * @struct ChainInfo
   * @brief Auxiliary data for tracking chain information for PAF output
   */
  struct ChainInfo {
    uint32_t chainId;    // Which chain this mapping belongs to
    uint16_t chainPos;   // Position within chain (1-based)
    uint16_t chainLen;   // Total mappings in this chain
  };
  
  typedef std::vector<ChainInfo> ChainInfoVector_t;
  
  /**
   * @struct MappingsWithChains
   * @brief Container for mappings with their associated chain information
   */
  struct MappingsWithChains {
    MappingResultsVector_t mappings;
    ChainInfoVector_t chainInfo;
  };

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
}

#endif
