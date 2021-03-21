/**
 * @file    computeAlignments.hpp
 * @brief   logic for generating alignments when given mashmap 
 *          mappings as input
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef COMPUTE_ALIGNMENTS_HPP 
#define COMPUTE_ALIGNMENTS_HPP

#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <zlib.h>
#include <cassert>
#include <thread>
#include <memory>

//Own includes
#include "align/include/align_types.hpp"
#include "align/include/align_parameters.hpp"
#include "map/include/base_types.hpp"
#include "map/include/commonFunc.hpp"

//External includes
#include "common/atomic_queue/atomic_queue.h"
#include "common/seqiter.hpp"
#include "common/progress.hpp"

#include "common/wflign/src/wflign_wfa.hpp"

namespace align
{

  long double float2phred(long double prob) {
    if (prob == 1)
        return 255;  // guards against "-0"
    long double p = -10 * (long double) log10(prob);
    if (p < 0 || p > 255) // int overflow guard
        return 255;
    else
        return p;
  }

  struct seq_record_t {
      MappingBoundaryRow currentRecord;
      std::string mappingRecordLine;
      std::shared_ptr<std::string> qSequence;
      seq_record_t(const MappingBoundaryRow& c, const std::string& r, const shared_ptr<std::string>& q)
          : currentRecord(c)
          , mappingRecordLine(r)
          , qSequence(q)
          { }
  };
  // load into this
  typedef atomic_queue::AtomicQueue<seq_record_t*, 2 << 16> seq_atomic_queue_t;
  // results into this, write out
  typedef atomic_queue::AtomicQueue<std::string*, 2 << 16> paf_atomic_queue_t;

  /**
   * @class     align::Aligner
   * @brief     compute alignments and generate sam output
   *            from mashmap mappings
   */
  class Aligner
  {
    private:

      //algorithm parameters
      const align::Parameters &param;

      refSequenceMap_t refSequences;

    public:

      /**
       * @brief                 constructor, also reads reference sequences
       * @param[in] p           algorithm parameters
       */
      Aligner(const align::Parameters &p) :
        param(p)
      {
        this->getRefSequences();
      }

      /**
       * @brief                 compute alignments
       */
      void compute()
      {
        this->computeAlignments();
      }

    private:

      /**
       * @brief                 parse and save all the reference sequences
       */
      void getRefSequences()
      {
        for(const auto &fileName : param.refSequences)
        {

#ifdef DEBUG
          std::cerr << "INFO, align::Aligner::getRefSequences, parsing reference sequences in file " << fileName << std::endl;
#endif

        seqiter::for_each_seq_in_file(
            fileName,
            [&](const std::string& seq_name,
                const std::string& seq) {
                // todo: offset_t is an 32-bit integer, which could cause problems
                skch::offset_t len = seq.length();

                skch::CommonFunc::makeUpperCase((char*)seq.c_str(), len);

                //seqId shouldn't already exist in our table
                assert(this->refSequences.count(seq_name) == 0);
                refSequences.emplace(seq_name, seq);
            });
        }
      }

      /**
       * @brief                 parse query sequences and mashmap mappings
       *                        to compute sequence alignments
       */
      void computeAlignments()
      {

          uint64_t total_seqs = 0;
          uint64_t total_alignment_length = 0;
          uint64_t total_paf_records = 0;
          for(const auto &fileName : param.querySequences) {
              std::ifstream mappingListStream(param.mashmapPafFile);
              std::string mappingRecordLine;
              MappingBoundaryRow currentRecord;
              seqiter::for_each_seq_in_file(
                  fileName,
                  [&](const std::string& qSeqId,
                      const std::string& _seq) {
                      ++total_seqs;
                      //total_seq_length += seq.size();
                      while(!mappingListStream.eof() && mappingRecordLine.empty()) {
                          std::getline(mappingListStream, mappingRecordLine);
                      }

                      if( !mappingRecordLine.empty() ) {
                          this->parseMashmapRow(mappingRecordLine, currentRecord);

                          if(currentRecord.qId == qSeqId) {
                              //auto q = new seq_record_t(currentRecord, mappingRecordLine, seq);
                              //seq_queue.push(q);
                              total_alignment_length += currentRecord.qEndPos - currentRecord.qStartPos;
                              ++total_paf_records;
                              //Check if more mappings have same query sequence id
                              while(std::getline(mappingListStream, mappingRecordLine)) {
                                  this->parseMashmapRow(mappingRecordLine, currentRecord);
                                  if(currentRecord.qId != qSeqId) {
                                      break;
                                  } else {
                                      total_alignment_length += currentRecord.qEndPos - currentRecord.qStartPos;
                                      ++total_paf_records;
                                  }
                              }
                          }
                      }
                  });
          }

          progress_meter::ProgressMeter progress(total_alignment_length, "[wfmash::align::computeAlignments] aligned");

          // input atomic queue
          seq_atomic_queue_t seq_queue;
          // output atomic queue
          paf_atomic_queue_t paf_queue;
          // flag when we're done reading
          std::atomic<bool> reader_done;
          reader_done.store(false);

          auto& nthreads = param.threads;

          // atomics to record if we're working or not
          std::vector<std::atomic<bool>> working(nthreads);
          for (auto& w : working) {
              w.store(true);
          }

          // reader picks up candidate alignments from input
          auto reader_thread =
              [&](void) {
                  //Parse query sequences
                  for(const auto &fileName : param.querySequences)
                  {
//#define DEBUG true
#ifdef DEBUG
                      std::cerr << "INFO, align::Aligner::computeAlignments, parsing query sequences in file " << fileName << std::endl;
#endif

                      //Open mashmap output file
                      std::ifstream mappingListStream(param.mashmapPafFile);
                      std::string mappingRecordLine;
                      MappingBoundaryRow currentRecord;

                      seqiter::for_each_seq_in_file(
                          fileName,
                          [&](const std::string& qSeqId,
                              const std::string& _seq) {
                              // copy our input into a shared ptr
                              std::shared_ptr<std::string> seq(new std::string(_seq));
                              // todo: offset_t is an 32-bit integer, which could cause problems
                              skch::offset_t len = seq->length();
                              // upper-case our input
                              skch::CommonFunc::makeUpperCase((char*)seq->c_str(), len);
                              // todo maybe this should change to some kind of unique pointer?
                              // something where we can GC it when we're done aligning to it
                              //std::string qSequence = seq;
                              //std::cerr << seq << std::endl;

                              //Check if all mapping records are processed already
                              while(!mappingListStream.eof() && mappingRecordLine.empty()) {
                                  //Read first record from mashmap output file during first iteration
                                  std::getline(mappingListStream, mappingRecordLine);
                              }

                              if( !mappingRecordLine.empty() ) {
                                  this->parseMashmapRow(mappingRecordLine, currentRecord);
                              
                                  //Check if mapping query id matches current query sequence id
                                  if(currentRecord.qId == qSeqId)
                                  {
                                      //Continue to read the next query sequence
                                      //continue;
                                      auto q = new seq_record_t(currentRecord, mappingRecordLine, seq);
                                      seq_queue.push(q);
                                  
                                      //Check if more mappings have same query sequence id
                                      while(std::getline(mappingListStream, mappingRecordLine))
                                      {
                                          this->parseMashmapRow(mappingRecordLine, currentRecord);
                                          
                                          if(currentRecord.qId != qSeqId)
                                          {
                                              //Break the inner loop to read query sequence
                                              break;
                                          }
                                          else
                                          {
                                              auto q = new seq_record_t(currentRecord, mappingRecordLine, seq);
                                              seq_queue.push(q);
                                          }
                                      }
                                  }
                              }
                          });
                          
                      mappingListStream.close();

                  }
                  reader_done.store(true);
              };

          // helper to check if we're still aligning
          auto still_working =
              [&](const std::vector<std::atomic<bool>>& working) {
                  bool ongoing = false;
                  for (auto& w : working) {
                      ongoing = ongoing || w.load();
                  }
                  return ongoing;
              };

          // writer, picks output from queue and writes it to our output stream
          std::ofstream outstrm(param.pafOutputFile);
          auto writer_thread =
              [&](void) {
                  while (true) {
                      std::string* paf_lines = nullptr;
                      if (!paf_queue.try_pop(paf_lines)
                          && !still_working(working)) {
                          break;
                      } else if (paf_lines != nullptr) {
                          outstrm << *paf_lines;
                          delete paf_lines;
                      } else {
                          std::this_thread::sleep_for(100ns);
                      }
                  }
              };

          // worker, takes candidate alignments and runs edlib alignment on them
          auto worker_thread = 
              [&](uint64_t tid,
                  std::atomic<bool>& is_working) {
                  is_working.store(true);
                  while (true) {
                      seq_record_t* rec = nullptr;
                      if (!seq_queue.try_pop(rec)
                          && reader_done.load()) {
                          break;
                      } else if (rec != nullptr) {
                          std::string* paf_rec
                              = new std::string(
                                  doAlignment(rec->currentRecord,
                                              rec->mappingRecordLine,
                                              rec->qSequence));
                          progress.increment(rec->currentRecord.qEndPos
                                             - rec->currentRecord.qStartPos);
                          if (paf_rec->size()) {
                              paf_queue.push(paf_rec);
                          } else {
                              delete paf_rec;
                          }
                          delete rec;
                      } else {
                          std::this_thread::sleep_for(100ns);
                      }
                  }
                  is_working.store(false);
              };

          // launch reader
          std::thread reader(reader_thread);
          // launch writer
          std::thread writer(writer_thread);
          // launch workers
          std::vector<std::thread> workers; workers.reserve(nthreads);
          for (uint64_t t = 0; t < nthreads; ++t) {
              workers.emplace_back(worker_thread,
                                   t,
                                   std::ref(working[t]));
          }

          // wait for reader and workers to complete
          reader.join();
          for (auto& worker : workers) {
              worker.join();
          }
          // and finally the writer
          writer.join();

          progress.finish();
          std::cerr << "[wfmash::align::computeAlignments] "
                    << "count of mapped reads = " << total_seqs
                    << ", total aligned bp = " << total_alignment_length << std::endl;

      }

      /**
       * @brief                         parse mashmap row sequence 
       * @param[in]   mappingRecordLine
       * @param[out]  currentRecord
       */
      inline void parseMashmapRow(const std::string &mappingRecordLine, MappingBoundaryRow &currentRecord)
      {
        std::stringstream ss(mappingRecordLine); // Insert the string into a stream
        std::string word; // Have a buffer string

        vector<std::string> tokens; // Create vector to hold our words

        while (ss >> word)
          tokens.push_back(word);

        //We expect and need at least these many values in a mashmap mapping
        assert(tokens.size() >= 9);

        //Save words into currentRecord
        {
          currentRecord.qId = tokens[0];
          currentRecord.qStartPos = std::stoi(tokens[2]);
          currentRecord.qEndPos = std::stoi(tokens[3]);
          currentRecord.strand = (tokens[4] == "+" ? skch::strnd::FWD : skch::strnd::REV);
          currentRecord.refId = tokens[5];
          currentRecord.rStartPos = std::stoi(tokens[7]);
          currentRecord.rEndPos = std::stoi(tokens[8]);
        }
      }

      /**
       * @brief                           compute alignment using edlib 
       * @param[in]   currentRecord       mashmap mapping parsed information
       * @param[in]   mappingRecordLine   mashmap mapping output raw string
       * @param[in]   qSequence           query sequence
       * @param[in]   outstrm             output stream
       */
      std::string doAlignment(MappingBoundaryRow &currentRecord,
                              const std::string &mappingRecordLine,
                              const std::shared_ptr<std::string> &qSequence) {

#ifdef DEBUG
        std::cerr << "INFO, align::Aligner::doAlignment, aligning mashmap record: " << mappingRecordLine << std::endl;
#endif

        //Define reference substring for this mapping
        const std::string &refId = currentRecord.refId;
        const char* refRegion = this->refSequences[refId].c_str();
        const auto& refSize = this->refSequences[refId].size();
        refRegion += currentRecord.rStartPos;
        skch::offset_t refLen = currentRecord.rEndPos - currentRecord.rStartPos;
        assert(refLen <= refSize);

        //Define query substring for this mapping
        const char* queryRegion = qSequence->c_str();  //initially point to beginning
        const auto& querySize = qSequence->size();
        skch::offset_t queryLen = currentRecord.qEndPos - currentRecord.qStartPos;
        queryRegion += currentRecord.qStartPos;

        char* queryRegionStrand = new char[queryLen+1];

        if(currentRecord.strand == skch::strnd::FWD) {
          strncpy(queryRegionStrand, queryRegion, queryLen);    //Copy the same string
        } else {
          skch::CommonFunc::reverseComplement(queryRegion, queryRegionStrand, queryLen); //Reverse complement
        }

        assert(queryLen <= querySize);

        //Compute alignment
        auto t0 = skch::Time::now();

#ifdef DEBUG
        std::cerr << "INFO, align::Aligner::doAlignment, WFA execution starting, query region length = " << queryLen
          << ", reference region length= " << refLen << ", edit distance limit= " << editDistanceLimit << std::endl; 
#endif

        std::stringstream output;
        // todo:
        // - toggle between wflign and regular alignment at some threshold (in wflign?)
        wflign::wavefront::wflign_affine_wavefront(
            output,
            true, // merge alignments
            !param.sam_format,
            currentRecord.qId, queryRegionStrand, querySize, currentRecord.qStartPos, queryLen,
            currentRecord.strand != skch::strnd::FWD,
            refId, refRegion, refSize, currentRecord.rStartPos, refLen,
            param.wflambda_segment_length,
            param.min_identity,
            param.wflambda_min_wavefront_length,
            param.wflambda_max_distance_threshold);

        delete [] queryRegionStrand;

        return output.str();
      }
  };
}


#endif
