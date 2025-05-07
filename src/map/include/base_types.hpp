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
    // IMPORTANT: This operator compares by coordinates, not by hash!
    bool operator <(const IntervalPoint& x) const {
      return std::tie(seqId, pos, side) 
        < std::tie(x.seqId, x.pos, x.side);
    }
  };

  // Type alias for the sorted vector of interval points
  using SortedIndexPoints = std::vector<IntervalPoint>;

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

  //Vector type for storing MappingResult objects (ownership)
  typedef std::vector<MappingResult> MappingResultsVector_t;

  //Vector type for views of MappingResult objects (non-owning pointers)
  typedef std::vector<MappingResult*> MappingResultsView_t;

  /**
   * @brief Create a view from a vector of MappingResult objects
   * @param mappings Vector of MappingResult objects (ownership)
   * @return Vector of pointers to MappingResult objects
   */
  inline MappingResultsView_t createViewFromMappings(MappingResultsVector_t& mappings) {
    MappingResultsView_t view;
    view.reserve(mappings.size());
    for (auto& mapping : mappings) {
      view.push_back(&mapping);
    }
    return view;
  }

  /**
   * @brief Add a new MappingResult to an owner vector and return a pointer to it
   * @param mappingsOwner Vector that will own the MappingResult
   * @param result The MappingResult to add
   * @return Pointer to the added MappingResult in the owner vector
   */
  inline MappingResult* addToOwner(MappingResultsVector_t& mappingsOwner, const MappingResult& result) {
    mappingsOwner.push_back(result);
    return &mappingsOwner.back();
  }

  /**
   * @brief Add a new MappingResult to an owner vector and add its pointer to a view
   * @param mappingsOwner Vector that will own the MappingResult
   * @param mappingsView View that will get a pointer to the new MappingResult
   * @param result The MappingResult to add
   */
  inline void addToOwnerAndView(MappingResultsVector_t& mappingsOwner, 
                               MappingResultsView_t& mappingsView, 
                               const MappingResult& result) {
    mappingsOwner.push_back(result);
    mappingsView.push_back(&mappingsOwner.back());
  }

  //Vector type for storing MinmerInfo
  typedef std::vector<MinmerInfo> MinVec_Type;

  //Helper functions for working with MappingResultsView_t
  
  /**
   * @brief Copy MappingResult objects from a view to an owner vector
   * @param view Source view containing pointers to MappingResult objects
   * @param owner Destination vector to own the copied MappingResult objects
   */
  inline void copyFromViewToOwner(const MappingResultsView_t& view, MappingResultsVector_t& owner) {
    owner.reserve(owner.size() + view.size());
    for (const auto* mapping_ptr : view) {
      owner.push_back(*mapping_ptr);
    }
  }

  /**
   * @brief Filter a view in-place using a predicate function
   * @param view Vector of pointers to MappingResult objects
   * @param predicate Function that takes a MappingResult pointer and returns bool
   */
  template<typename Predicate>
  inline void filterViewInPlace(MappingResultsView_t& view, Predicate predicate) {
    view.erase(
      std::remove_if(view.begin(), view.end(), 
                    [&predicate](const MappingResult* m) { return !predicate(m); }),
      view.end()
    );
  }

  /**
   * @brief Mark mappings in a view as discarded
   * @param view Vector of pointers to MappingResult objects to be discarded
   */
  inline void markViewAsDiscarded(const MappingResultsView_t& view) {
    for (auto* mapping_ptr : view) {
      mapping_ptr->discard = 1;
    }
  }

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
