/**
 * @file    base_types.hpp
 * @brief   Critical type definitions for mapping algorithm
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef BASE_TYPES_MAP_HPP
#define BASE_TYPES_MAP_HPP

#include <tuple>
#include <vector>

namespace skch
{
  typedef uint32_t hash_t;    //hash type
  typedef int offset_t;       //position within sequence
  typedef int seqno_t;        //sequence counter in file
  typedef int32_t strand_t;   //sequence strand

  //C++ timer
  typedef std::chrono::high_resolution_clock Time;

  //Information about each minimizer
  struct MinimizerInfo
  {
    hash_t hash;                              //hash value
    seqno_t seqId;                            //sequence or contig id
    offset_t wpos;                            //First (left-most) window position when the minimizer is saved
    strand_t strand;                          //strand information
    //double order;

    //Lexographical less than comparison
    bool operator <(const MinimizerInfo& x) {
      return std::tie(hash, seqId, wpos, strand)
        < std::tie(x.hash, x.seqId, x.wpos, x.strand);
    }

    //Lexographical equality comparison
    bool operator ==(const MinimizerInfo& x) {
      return std::tie(hash, seqId, wpos, strand)
        == std::tie(x.hash, x.seqId, x.wpos, x.strand);
    }

    bool operator !=(const MinimizerInfo& x) {
      return std::tie(hash, seqId, wpos, strand)
        != std::tie(x.hash, x.seqId, x.wpos, x.strand);
    }

    static bool equalityByHash(const MinimizerInfo& x, const MinimizerInfo& y) {
      return x.hash == y.hash;
    }

    static bool lessByHash(const MinimizerInfo& x, const MinimizerInfo& y) {
      return x.hash < y.hash;
    }

  };

  //Type for map value type used for
  //L1 stage lookup index
  struct MinimizerMetaData
  {
    seqno_t seqId;          //sequence or contig id
    offset_t wpos;          //window position (left-most window)
    strand_t strand;        //strand information

    bool operator <(const MinimizerMetaData& x) const {
      return std::tie(seqId, wpos, strand)
        < std::tie(x.seqId, x.wpos, x.strand);
    }
  };

  typedef hash_t MinimizerMapKeyType;
  typedef std::vector<MinimizerMetaData> MinimizerMapValueType;

  //Metadata recording for contigs in the reference DB
  struct ContigInfo
  {
    std::string name;       //Name of the sequence
    offset_t len;           //Length of the sequence
  };

  //Label tags for strand information
  enum strnd : strand_t
  {
    FWD = 1,
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
    int blockLength;                                    //the block length of the mapping
    float nucIdentity;                                  //calculated identity
    float nucIdentityUpperBound;                        //upper bound on identity (90% C.I.)
    int sketchSize;                                     //sketch size
    int conservedSketches;                              //count of conserved sketches
    strand_t strand;                                    //strand
    int approxMatches;                                  //the approximate number of matches in the alignment

                                                        //--for split read mapping

    int n_merged;                                       // how many mappings we've merged into this one
    offset_t splitMappingId;                            // To identify split mappings that are chained
    uint8_t discard;                                    // set to 1 for deletion
    bool selfMapFilter;                                 // set to true if a long-to-short mapping in all-vs-all mode (we report short as the query)

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

  //Container to save copy of kseq object
  struct InputSeqContainer
  {
    seqno_t seqCounter;                 //sequence counter
    offset_t len;                       //sequence length
    std::string seq;                    //sequence string
    std::string seqName;                //sequence id

    /*
     * @brief               constructor
     * @param[in] kseq_seq  complete read or reference sequence
     * @param[in] kseq_id   sequence id name
     * @param[in] len       length of sequence
     */
      InputSeqContainer(const std::string& s, const std::string& id, seqno_t seqcount)
          : seq(s)
          , seqName(id)
          , len(s.length())
          , seqCounter(seqcount) { }
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
  template <typename MinimizerVec>
    struct QueryMetaData
    {
      char *seq;                          //query sequence pointer
      seqno_t seqCounter;                 //query sequence counter
      offset_t len;                       //length of this query sequence
      offset_t fullLen;                   //length of the full sequence it derives from
      int sketchSize;                     //sketch size
      std::string seqName;                //sequence name
      MinimizerVec minimizerTableQuery;   //Vector of minimizers in the query
    };
}

#endif
