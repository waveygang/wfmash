/**
 * @file    filter.hpp
 * @brief   implements the routines to filter mappings
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef FILTER_MAP_HPP
#define FILTER_MAP_HPP

#include <vector>
#include <algorithm>
#include <set>
#include <fstream>
#include <zlib.h>
#include <iostream>
#include <sstream>

//Own includes
#include "map/include/base_types.hpp"
#include "map/include/map_parameters.hpp"
#include "common/progress.hpp"

//External includes

namespace skch
{
  /**
   * @namespace skch::CommonFunc
   * @brief     Implements routines to filter mappings
   */
  namespace Filter
  {
    /**
     * @namespace skch::filter::query
     * @brief     filter routines (best for query sequence)
     */

    namespace query
    {
      //helper functions for executing plane sweep over query sequence
      struct Helper
      {
        MappingResultsVector_t &vec;

        Helper(MappingResultsVector_t &v) : vec(v) {}

        double get_score(const int x) const {
            if (vec[x].blockLength <= 0 || vec[x].blockNucIdentity() <= 0) {
                return std::numeric_limits<double>::lowest();
            }
            return vec[x].blockNucIdentity() * std::log(static_cast<double>(vec[x].blockLength));
        }

        //Greater than comparison by score and begin position
        //used to define order in BST
        bool operator ()(const int x, const int y) const {

          assert(x < vec.size());
          assert(y < vec.size());

          auto x_score = get_score(x);
          auto y_score = get_score(y);

          return std::tie(x_score, vec[x].queryStartPos, vec[x].refSeqId) > std::tie(y_score, vec[y].queryStartPos, vec[y].refSeqId);
        }

        //Greater than comparison by score
        bool greater_score(const int x, const int y) const {

          assert(x < vec.size());
          assert(y < vec.size());

          auto x_score = get_score(x);
          auto y_score = get_score(y);

          return x_score > y_score;
        }

        // compute the overlap of the two mappings
        double get_overlap(const int x, const int y) const {
            offset_t overlap_start = std::max(vec[x].queryStartPos, vec[y].queryStartPos);
            offset_t overlap_end = std::min(vec[x].queryEndPos(), vec[y].queryEndPos());
            offset_t overlap_length = std::max(0, static_cast<int>(overlap_end - overlap_start));
            offset_t x_length = vec[x].queryEndPos() - vec[x].queryStartPos;
            offset_t y_length = vec[y].queryEndPos() - vec[y].queryStartPos;
            return static_cast<double>(overlap_length) / std::min(x_length, y_length);
        }

        /*
         * @brief                         mark the mappings with maximum score as good (on query seq)
         * @tparam          Type          std::set type to save mappings (sweep line status container)
         * @param[in/out]   L             container with mappings
         */
        template <typename Type>
        inline void markGood(Type &L, int secondaryToKeep, bool dropRand, double overlapThreshold)
          {
            //first segment in the set order
            auto beg = L.begin();

            // count how many secondary alignments we keep
            int kept = 0;

            auto it = L.begin();
            for( ; it != L.end(); it++)
            {
                if ((this->greater_score(*beg, *it) || !vec[*it].discard()) && kept > secondaryToKeep) {
                    break;
                }

                vec[*it].setDiscard(false);
                ++kept;
            }
            auto kit = it;

            // Skip overlap checking if threshold is 1.0 (allow all overlaps)
            if (overlapThreshold < 1.0) {
                // Check for overlaps and mark bad if necessary
                for ( ; it != L.end(); it++) {
                    if (it == L.begin()) continue;
                    int idx = *it;
                    for (auto it2 = L.begin(); it2 != kit; it2++) {
                        double overlap = get_overlap(idx, *it2);
                        if (overlap > overlapThreshold) {
                            vec[idx].setOverlapped(true);  // Mark as bad if overlaps more than threshold
                            vec[idx].setDiscard(true);
                            break;
                        }
                    }
                }
            }

            // check for the case where there are multiple best mappings > secondaryToKeep
            // which have the same score
            // we will hash the mapping struct and keep the one with the secondaryToKeep with the lowest hash value
            if (kept > secondaryToKeep && dropRand) 
            {
              // we will use hashes of the mapping structs to break ties
              // first we'll make a vector of the mappings including the hashes
              std::vector<std::tuple<double, size_t, MappingResult*>> score_and_hash; // The tuple is (score, hash, pointer to the mapping)
              for(auto it = L.begin(); it != L.end(); it++)
              {
                  if(!vec[*it].discard())
                  {
                      score_and_hash.emplace_back(get_score(*it), vec[*it].hash(), &vec[*it]);
                  }
              }
              // now we'll sort the vector by score and hash
              std::sort(score_and_hash.begin(), score_and_hash.end(), std::greater{});
              // reset kept counter
              kept = 0;
              for (auto& x : score_and_hash) {
                  std::get<2>(x)->setDiscard(true);
              }
              // now we mark the best to keep
              for (auto& x : score_and_hash) {
                  if (kept > secondaryToKeep) {
                      break;
                  }
                  std::get<2>(x)->setDiscard(false);
                  ++kept;
              }
            }
          }
      };

     /**
       * @brief                       filter mappings (best for query sequence)
       * @details                     evaluate best cover mapping for each base pair
       * @param[in/out] readMappings  Mappings computed by Mashmap
       */
      template <typename VecIn>
      void liFilterAlgorithm(VecIn &readMappings, int secondaryToKeep, bool dropRand, double overlapThreshold, progress_meter::ProgressMeter& progress)
        {
          if(readMappings.size() <= 1)
            return;

          //Initially mark all mappings as bad
          //Maintain the order of this vector till end of this function
          std::for_each(readMappings.begin(), readMappings.end(), [&](MappingResult &e){ 
            e.setDiscard(true); 
            e.setOverlapped(false); 
          });

          //Initialize object of Helper struct
          Helper obj (readMappings);

          //Plane sweep status
          //binary search tree of segment ids, ordered by their scores
          std::set <int, Helper> bst (obj);

          //Event point schedule
          //vector of triplets <position, event type, segment id>
          typedef std::tuple<offset_t, int, int> eventRecord_t;
          std::vector <eventRecord_t>  eventSchedule (2*readMappings.size());

          for(int i = 0; i < readMappings.size(); i++)
          {
            eventSchedule.emplace_back (readMappings[i].queryStartPos, event::BEGIN, i);
            eventSchedule.emplace_back (readMappings[i].queryEndPos(), event::END, i);
          }

          std::sort(eventSchedule.begin(), eventSchedule.end());

          //Execute the plane sweep algorithm
          for(auto it = eventSchedule.begin(); it!= eventSchedule.end();)
          {
            //Find events that correspond to current position
            auto it2 = std::find_if(it, eventSchedule.end(), [&](const eventRecord_t &e)
                                    {
                                      return std::get<0>(e) != std::get<0>(*it);
                                    });

            //update sweep line status by adding/removing segments
            std::for_each(it, it2, [&](const eventRecord_t &e)
                                    {
                                      int idx = std::get<2>(e);
                                      if (std::get<1>(e) == event::BEGIN)
                                        bst.insert (idx);
                                      else
                                        bst.erase (idx);
                                    });

            //mark mappings as good
            obj.markGood(bst, secondaryToKeep, dropRand, overlapThreshold);

            it = it2;
          }

          //Remove bad mappings
          readMappings.erase(
              std::remove_if(readMappings.begin(), readMappings.end(), [&](MappingResult &e){ 
                return e.discard() || e.overlapped();
              }),
              readMappings.end());
        }

      /**
       * @brief                       filter mappings (best for query sequence)
       * @details                     evaluate best N unmerged mappings for each position, assumes non-overlapping mappings in query
       * @param[in/out] readMappings  Mappings computed by Mashmap
       */
      template <typename VecIn>
      void indexedFilterAlgorithm(VecIn &readMappings, int secondaryToKeep)
        {
          if(readMappings.size() <= 1)
            return;

          //Initially mark all mappings as bad
          //Maintain the order of this vector till end of this function
          std::for_each(readMappings.begin(), readMappings.end(), [&](MappingResult &e){ e.setDiscard(true); });

          //Initialize object of Helper struct
          Helper obj (readMappings);

          //Event point schedule
          //vector of triplets <position, event type, segment id>
          typedef std::tuple<offset_t, double, int, int> eventRecord_t;
          std::vector <eventRecord_t>  eventSchedule (2*readMappings.size());

          for(int i = 0; i < readMappings.size(); i++) {
              eventSchedule.emplace_back (readMappings[i].queryStartPos, obj.get_score(i), event::BEGIN, i);
              eventSchedule.emplace_back (readMappings[i].queryEndPos, 0, event::END, i); // end should not be preferred
          }

          std::sort(eventSchedule.begin(), eventSchedule.end());

          //Execute the plane sweep algorithm
          for(auto it = eventSchedule.begin(); it!= eventSchedule.end();)
          {
            //Find events that correspond to current position
            auto it2 = std::find_if(it, eventSchedule.end(), [&](const eventRecord_t &e)
                                    {
                                      return std::get<0>(e) != std::get<0>(*it);
                                    });

            //mark best secondaryToKeep+1 mappings as good
            int kept = 0;
            std::for_each(it, it2, [&](const eventRecord_t &e)
                                    {
                                        if (std::get<2>(e) == event::BEGIN && kept <= secondaryToKeep) {
                                            obj.vec[std::get<3>(e)].setDiscard(false);
                                            ++kept;
                                        }
                                    });

            it = it2;
          }

          //Remove bad mappings
          readMappings.erase(
              std::remove_if(readMappings.begin(), readMappings.end(), [&](MappingResult &e){ return e.discard(); }),
              readMappings.end());
        }

      /**
       * @brief                          filter mappings (best for query sequence)
       * @param[in/out] readMappings     Mappings computed by Mashmap (post merge step)
       * @param[in]     secondaryToKeep  How many mappings in addition to the best to keep
       * @param[in]     dropRand         If multiple mappings have the same score, drop randomly
       *                                 until we only have secondaryToKeep secondary mappings
       */
      template <typename VecIn>
      void filterMappings(VecIn &readMappings, int secondaryToKeep, bool dropRand, double overlapThreshold, progress_meter::ProgressMeter& progress)
      {
          //Apply the main filtering algorithm to ensure the best mappings across complete axis
          liFilterAlgorithm(readMappings, secondaryToKeep, dropRand, overlapThreshold, progress);
      }

     /**
       * @brief                       filter mappings (best for query sequence)
       * @param[in/out] readMappings  Mappings computed by Mashmap (post merge step)
       */
      template <typename VecIn>
      void filterUnmergedMappings(VecIn &readMappings, int secondaryToKeep, progress_meter::ProgressMeter& progress)
      {
          //Apply a simple filtering algorithm that keeps the best secondaryToKeep+1 mappings per position
          indexedFilterAlgorithm(readMappings, secondaryToKeep, progress);
      }
    } //End of query namespace

    namespace ref
    {

      //Event point schedule
      //vector of triplets <ref sequence id, ref seq offset, event type, segment id>
      typedef std::tuple <seqno_t, offset_t, int, int> eventRecord_t;

      //helper functions for executing plane sweep over reference sequence
      struct Helper
      {
        MappingResultsVector_t &vec;

        Helper(MappingResultsVector_t &v) : vec(v) {}

        double get_score(const int x) const {return vec[x].blockNucIdentity() * log(vec[x].blockLength); }

        //Greater than comparison by score and begin position
        //used to define order in BST
        bool operator ()(const int x, const int y) const {

          assert(x < vec.size());
          assert(y < vec.size());

          auto x_score = get_score(x);
          auto y_score = get_score(y);

          return std::tie(x_score, vec[x].refStartPos) > std::tie(y_score, vec[y].refStartPos);
        }

        //Greater than comparison by score
        bool greater_score(const int x, const int y) const {

          assert(x < vec.size());
          assert(y < vec.size());

          auto x_score = get_score(x);
          auto y_score = get_score(y);

          return x_score > y_score;
        }

        // compute the overlap of the two mappings
        double get_overlap(const int x, const int y) const {
            offset_t overlap_start = std::max(vec[x].refStartPos, vec[y].refStartPos);
            offset_t overlap_end = std::min(vec[x].refEndPos(), vec[y].refEndPos());
            offset_t overlap_length = std::max(0, static_cast<int>(overlap_end - overlap_start));
            offset_t x_length = vec[x].refEndPos() - vec[x].refStartPos;
            offset_t y_length = vec[y].refEndPos() - vec[y].refStartPos;
            return static_cast<double>(overlap_length) / std::min(x_length, y_length);
        }

        /**
         * @brief                         mark the mappings with maximum score as good (on query seq)
         * @tparam          Type          std::set type to save mappings (sweep line status container)
         * @param[in/out]   L             container with mappings
         */
          template <typename Type>
          inline void markGood(Type &L, int secondaryToKeep, bool dropRand, double overlapThreshold)
          {
            //first segment in the set order
            auto beg = L.begin();

            // count how many secondary alignments we keep
            int kept = 0;

            auto it = L.begin();
            for( ; it != L.end(); it++)
            {
                if ((this->greater_score(*beg, *it) || !vec[*it].discard()) && kept > secondaryToKeep) {
                    break;
                }

                vec[*it].setDiscard(false);
                ++kept;
            }
            auto kit = it;

            // Skip overlap checking if threshold is 1.0 (allow all overlaps)
            if (overlapThreshold < 1.0) {
                // Check for overlaps and mark bad if necessary
                for ( ; it != L.end(); it++) {
                    if (it == L.begin()) continue;
                    int idx = *it;
                    for (auto it2 = L.begin(); it2 != kit; it2++) {
                        if (get_overlap(idx, *it2) > overlapThreshold) {
                            vec[idx].setOverlapped(true);  // Mark as bad if overlaps more than threshold
                            vec[idx].setDiscard(true);
                            break;
                        }
                    }
                }
            }

            // check for the case where there are multiple best mappings > secondaryToKeep
            // which have the same score
            // we will hash the mapping struct and keep the one with the secondaryToKeep with the lowest hash value
            if (kept > secondaryToKeep && dropRand) 
            {
              // we will use hashes of the mapping structs to break ties
              // first we'll make a vector of the mappings including the hashes
              std::vector<std::tuple<double, size_t, MappingResult*>> score_and_hash; // The tuple is (score, hash, pointer to the mapping)
              for(auto it = L.begin(); it != L.end(); it++)
              {
                  if(!vec[*it].discard())
                  {
                      score_and_hash.emplace_back(get_score(*it), vec[*it].hash(), &vec[*it]);
                  }
              }
              // now we'll sort the vector by score and hash
              std::sort(score_and_hash.begin(), score_and_hash.end(), std::greater{});
              // reset kept counter
              kept = 0;
              for (auto& x : score_and_hash) {
                  std::get<2>(x)->setDiscard(true);
              }
              // now we mark the best to keep
              for (auto& x : score_and_hash) {
                  if (kept > secondaryToKeep) {
                      break;
                  }
                  std::get<2>(x)->setDiscard(false);
                  ++kept;
              }
            }
          }

        /**
         * @brief                       advance position by one over reference sequence(s)
         * @param[in/out] eventRecord   event record containing end point of segment
         * @param[in]     refsketch     reference index class object
         */
        void refPosDoPlusOne(eventRecord_t &eventRecord, const skch::SequenceIdManager &idManager)
        {
          seqno_t currentSeqId = std::get<0>(eventRecord);
          offset_t currentSeqOffSet = std::get<1>(eventRecord);

          //if offset is at the end of reference sequence, shift to next
          if(currentSeqOffSet == idManager.getSequenceLength(currentSeqId) - 1)
          {
            std::get<0>(eventRecord) += 1;    //shift id by 1
            std::get<1>(eventRecord) = 0;
          }
          else
            std::get<1>(eventRecord) += 1;    //shift offset by 1:
        }
      };

      /**
       * @brief                       filter mappings (best for reference sequence)
       * @param[in/out] readMappings  Mappings computed by Mashmap (post merge step)
       * @param[in]     refsketch     reference index class object, used to determine ref sequence lengths
       */
      template <typename VecIn>
      void filterMappings(VecIn &readMappings, const skch::SequenceIdManager &idManager, uint16_t secondaryToKeep, bool dropRand, double overlapThreshold)
        {
          if(readMappings.size() <= 1)
            return;

          //Initially mark all mappings as bad
          //Maintain the order of this vector till end of this function
          std::for_each(readMappings.begin(), readMappings.end(), [&](MappingResult &e){ e.setDiscard(true); });

          //Initialize object of Helper struct
          Helper obj (readMappings);

          //Plane sweep status
          //binary search tree of segment ids, ordered by their scores
          std::set <int, Helper> bst (obj);

          //Event point schedule
          //vector of triplets <position, event type, segment id>
          std::vector <eventRecord_t>  eventSchedule (2*readMappings.size());

          for(int i = 0; i < readMappings.size(); i++)
          {
            eventSchedule.emplace_back (readMappings[i].refSeqId, readMappings[i].refStartPos, event::BEGIN, i);

            eventRecord_t endEvent = std::make_tuple(readMappings[i].refSeqId, readMappings[i].refEndPos(), event::END, i);
            //add one to above coordinate
            obj.refPosDoPlusOne(endEvent, idManager);
            eventSchedule.push_back (endEvent);
          }

          std::sort(eventSchedule.begin(), eventSchedule.end());

          //Execute the plane sweep algorithm
          for(auto it = eventSchedule.begin(); it!= eventSchedule.end();)
          {
            //Find events that correspond to current position
            auto it2 = std::find_if(it, eventSchedule.end(), [&](const eventRecord_t &e)
                                    {
                                      return std::tie(std::get<0>(e), std::get<1>(e)) != std::tie(std::get<0>(*it), std::get<1>(*it));
                                    });

            //update sweep line status by adding/removing segments
            std::for_each(it, it2, [&](const eventRecord_t &e)
                                    {
                                      if (std::get<2>(e) == event::BEGIN)
                                        bst.insert (std::get<3>(e));
                                      else
                                        bst.erase (std::get<3>(e));
                                    });

            //mark mappings as good
            obj.markGood(bst, secondaryToKeep, dropRand, overlapThreshold);

            it = it2;
          }

          //Remove bad mappings
          readMappings.erase(
              std::remove_if(readMappings.begin(), readMappings.end(), [&](MappingResult &e){ return e.discard(); }),
              readMappings.end());
        }
    } //End of reference namespace
  }
}

#endif
