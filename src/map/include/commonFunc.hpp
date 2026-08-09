/**
 * @file    commonFunc.hpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef COMMON_FUNC_HPP
#define COMMON_FUNC_HPP

#include <vector>
#include <map>
#include <array>
#include <algorithm>
#include <deque>
#include <cmath>
#include <fstream>
#include <limits>
#include <numeric>
#include <queue>
#include <sstream>

//Own includes
#include "map/include/map_parameters.hpp"

//External includes
#include "common/murmur3.h"
#include "common/prettyprint.hpp"
#include "common/ankerl/unordered_dense.hpp"

//#include "assert.hpp"

namespace skch {
    /**
     * @namespace skch::CommonFunc
     * @brief     Implements frequently used common functions
     */
    namespace CommonFunc {
        //seed for murmerhash
        const int seed = 42;

        // Pivot to keep track of sketch border
        template <typename I>
        struct Pivot {
            I p;
            int64_t rank;
        };

        /**
         * @brief   reverse complement of kmer (borrowed from mash)
         * @note    assumes dest is pre-allocated
         */
        // Complement lookup table: A<->T, C<->G, everything else unchanged
        // (byte-identical to the per-base switch, but branchless).
        inline const std::array<char, 256> complement_table = [] {
            std::array<char, 256> t{};
            for (int i = 0; i < 256; ++i) t[i] = (char)i;
            t[(unsigned char)'A'] = 'T'; t[(unsigned char)'C'] = 'G';
            t[(unsigned char)'G'] = 'C'; t[(unsigned char)'T'] = 'A';
            return t;
        }();

        inline void reverseComplement(const char *src, char *dest, int length) {
            for (int i = 0; i < length; i++) {
                dest[length - i - 1] = complement_table[(unsigned char)src[i]];
            }
        }
        // Crazy hack char table to test for canonical bases
    constexpr int valid_dna[127] = {
        1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 0, 1, 0, 1, 1, 1,
        0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 0, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 0, 1, 0, 1,
        1, 1, 0, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 0, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1
    };

    /**
     * @brief               convert DNA or AA alphabets to upper case, converting non-canonical DNA bases to N
     * @param[in]   seq     pointer to input sequence
     * @param[in]   len     length of input sequence
     */
        // Precomputed uppercase+validate table: uppercase a-z, then map any non-ACGT
        // to 'N'. Byte-identical to the branch version for ASCII input; also defines
        // bytes >=127 (previously an out-of-bounds read of valid_dna) as 'N'.
        inline const std::array<char, 256> dna_clean_table = [] {
            std::array<char, 256> t{};
            for (int c = 0; c < 256; ++c) {
                int u = (c > 96 && c < 123) ? c - 32 : c;   // uppercase a-z
                bool invalid = (u >= 0 && u < 127) ? (bool)valid_dna[u] : true;
                t[c] = invalid ? 'N' : (char)u;
            }
            return t;
        }();

        inline void makeUpperCaseAndValidDNA(char *seq, offset_t len) {
            for (offset_t i = 0; i < len; i++) {
                seq[i] = dna_clean_table[(unsigned char)seq[i]];
            }
        }

//        /**
//       * @brief               convert DNA or AA alphabets to upper case
//       * @param[in]   seq     pointer to input sequence
//       * @param[in]   len     length of input sequence
//       */
//        inline void makeUpperCase(char *seq, offset_t len) {
//            for (int i = 0; i < len; i++) {
//                if (seq[i] > 96 && seq[i] < 123) {
//                    seq[i] -= 32;
//                }
//            }
//        }
//
//        /**
//         * @brief               convert non-canonical DNA bases to N
//         * @param[in]   seq     pointer to input sequence
//         * @param[in]   len     length of input sequence
//         */
//        inline void makeValidDNA(char *seq, offset_t len) {
//            for (int i = 0; i < len; i++) {
//                if (rkmh::valid_dna[seq[i]]) {
//                    seq[i] = 'N';
//                }
//            }
//        }

        /**
         * @brief   hashing kmer string (borrowed from mash)
         */
        inline hash_t getHash(const char *seq, int length) {
            char data[16];
            MurmurHash3_x64_128(seq, length, seed, data);

            hash_t hash;

            hash = *((hash_t *) data);

            return hash;
        }

        /**
         * @brief		takes hash value of kmer and adjusts it based on kmer's weight
         *					this value will determine its order for minimizer selection
         * @details	this is inspired from Chum et al.'s min-Hash and tf-idf weighting
         */
//        static inline double applyWeight(char* kmer, int kmer_size, hash_t kmer_hash, const std::unordered_set<std::string>& high_freq_kmers) {
//            double x = kmer_hash * 1.0 / UINT32_MAX;  //bring it within [0, 1]
//            //assert (x >= 0.0 && x <= 1.0);
//
//            std::string kmer_str(kmer, kmer_size);
//            if (high_freq_kmers.count(kmer_str) > 0) {
//                /* downweigting by a factor of 8 */
//                /* further aggressive downweigting may affect accuracy */
//                double p2 = x*x;
//                double p4 = p2 * p2;
//                return - 1.0 * (p4 * p4);
//            }
//            return -1.0 * x;
//
//            //range of returned value is between [-1,0]
//            //we avoid adding one for better double precision
//        }


        /**
         * @brief       Compute the minimum s kmers for a string.
         * @param[out]  minmerIndex     container storing sketched Kmers 
         * @param[in]   seq                 pointer to input sequence
         * @param[in]   len                 length of input sequence
         * @param[in]   kmerSize
         * @param[in]   s                   sketch size. 
         * @param[in]   seqCounter          current sequence number, used while saving the position of minimizer
         */
        template <typename T>
          inline void sketchSequence(
              std::vector<T> &minmerIndex, 
              char* seq, 
              offset_t len,
              int kmerSize, 
              int alphabetSize,
              int sketchSize,
              seqno_t seqCounter)
        {
          makeUpperCaseAndValidDNA(seq, len);

          //Compute reverse complement of seq (buffer reused across fragments)
          thread_local std::vector<char> seqRevBuf;
          seqRevBuf.resize(len);
          char* seqRev = seqRevBuf.data();

          if(alphabetSize == 4) //not protein
            CommonFunc::reverseComplement(seq, seqRev, len);

          // TODO cleanup
          thread_local ankerl::unordered_dense::map<hash_t, MinmerInfo> sketched_vals;
          sketched_vals.clear();
          sketched_vals.reserve(sketchSize + 1);
          thread_local std::vector<hash_t> sketched_heap;
          sketched_heap.clear();
          sketched_heap.reserve(sketchSize+1);
            
          // Get distance until last "N"
          int ambig_kmer_count = 0;
          for (int i = kmerSize - 1; i >= 0; i--)
          {
            if (seq[i] == 'N')
            {
                ambig_kmer_count = i+1;
                break;
            }    
          } 

          for(offset_t i = 0; i < len - kmerSize + 1; i++)
          {

            if (seq[i+kmerSize-1] == 'N')
            {
              ambig_kmer_count = kmerSize;
            }
            //Hash kmers
            hash_t hashFwd = CommonFunc::getHash(seq + i, kmerSize); 
            hash_t hashBwd;

            if(alphabetSize == 4)
              hashBwd = CommonFunc::getHash(seqRev + len - i - kmerSize, kmerSize);
            else  //proteins
              hashBwd = std::numeric_limits<hash_t>::max();   //Pick a dummy high value so that it is ignored later

            //Consider non-symmetric kmers only
            if(hashBwd != hashFwd && ambig_kmer_count == 0)
            {
              //Take minimum value of kmer and its reverse complement
              hash_t currentKmer = std::min(hashFwd, hashBwd);

              //Check the strand of this minimizer hash value
              auto currentStrand = hashFwd < hashBwd ? strnd::FWD : strnd::REV;

              if (sketched_heap.size() < sketchSize || currentKmer <= sketched_heap.front())
              {
                if (sketched_heap.empty() || sketched_vals.find(currentKmer) == sketched_vals.end()) 
                {

                  // Add current hash to heap
                  if (sketched_vals.size() < sketchSize || currentKmer < sketched_heap.front())  
                  {
                      sketched_vals[currentKmer] = MinmerInfo{currentKmer, i, i, seqCounter, currentStrand};
                      sketched_heap.push_back(currentKmer);
                      std::push_heap(sketched_heap.begin(), sketched_heap.end());
                  }

                  // Remove one if too large
                  if (sketched_vals.size() > sketchSize) 
                  {
                      sketched_vals.erase(sketched_heap[0]);
                      std::pop_heap(sketched_heap.begin(), sketched_heap.end());
                      sketched_heap.pop_back();
                  }
                } 
                else 
                {
                  // TODO these sketched values might never be useful, might save memory by deleting
                  // extend the length of the window
                  auto& sketchedVal = sketched_vals[currentKmer];
                  sketchedVal.wpos_end = i;
                  sketchedVal.strand += currentStrand == strnd::FWD ? 1 : -1;
                }
              }
            }
            if (ambig_kmer_count > 0)
            {
              ambig_kmer_count--;
            }
          }

          minmerIndex.resize(sketched_heap.size());
          for (auto rev_it = minmerIndex.rbegin(); rev_it != minmerIndex.rend(); rev_it++)
          {
            *rev_it = (std::move(sketched_vals[sketched_heap.front()]));
            (*rev_it).strand = (*rev_it).strand > 0 ? strnd::FWD : ((*rev_it).strand == 0 ? strnd::AMBIG : strnd::REV);

            std::pop_heap(sketched_heap.begin(), sketched_heap.end());
            sketched_heap.pop_back();
          }
          return;
        }
        

        /**
         * @brief       Compute winnowed minmers from a given sequence and add to the index
         * @param[out]  minmerIndex  table storing minmers and their position as we compute them
         * @param[in]   seq             pointer to input sequence
         * @param[in]   len             length of input sequence
         * @param[in]   kmerSize
         * @param[in]   windowSize
         * @param[in]   sketchSize      sketch size. 
         * @param[in]   seqCounter      current sequence number, used while saving the position of minimizer
         */
        template <typename T>
          inline void addMinmers(std::vector<T> &minmerIndex, 
              char* seq, offset_t len,
              int kmerSize, 
              int windowSize,
              int alphabetSize,
              int sketchSize,
              seqno_t seqCounter)
          {
            /**
             * Double-ended queue (saves minimum at front end)
             * Saves pair of the minimizer and the position of hashed kmer in the sequence
             * Position of kmer is required to discard kmers that fall out of current window
             */
            std::deque< std::tuple<hash_t, strand_t, offset_t> > Q;

            // One sliding-window entry per distinct in-window hash: the open
            // minmer plus a FIFO of its in-window occurrences (vector + head
            // index instead of std::deque; live range is [occHead, occ.size())).
            struct WindowEntry {
              MinmerInfo mi;
              std::vector<KmerInfo> occ;
              uint32_t occHead = 0;
              size_t count() const { return occ.size() - occHead; }
            };

            // Sort by hash, then by position
            constexpr auto KIHeap_cmp = [](KmerInfo& a, KmerInfo& b)
              {return std::tie(a.hash, a.pos) > std::tie(b.hash, b.pos);};
            // The former std::map<hash, (MinmerInfo, deque)> split into an O(1)
            // hash map plus an ascending vector of the in-window hashes: max is
            // back(), the final sweep scans in order, insert/erase are 8-byte
            // memmoves instead of tree rebalances.
            thread_local ankerl::unordered_dense::map<hash_t, WindowEntry> windowMap;
            windowMap.clear();
            thread_local std::vector<hash_t> sortedHashes;
            sortedHashes.clear();
            thread_local std::vector<KmerInfo> heapWindow;
            heapWindow.clear();

            makeUpperCaseAndValidDNA(seq, len);

            //Compute reverse complement of the whole sequence once. Hashing the
            //corresponding RC slice below is byte-identical to per-kmer RC but avoids
            //recomputing an O(kmerSize) reverse-complement at every position.
            std::unique_ptr<char[]> seqRev(new char[len]);
            if(alphabetSize == 4) //not protein
              CommonFunc::reverseComplement(seq, seqRev.get(), len);

            // Get distance until last "N"
            int ambig_kmer_count = 0;


            for(offset_t i = 0; i < len - kmerSize + 1; i++)
            {
              //The serial number of current sliding window
              //First valid window appears when i = windowSize - 1
              offset_t currentWindowId = i + kmerSize - windowSize;

              // Remove expired kmers from heap
              if (heapWindow.size() > 2*windowSize)
              {
                heapWindow.erase(
                    std::remove_if(
                      heapWindow.begin(), 
                      heapWindow.end(),
                      [currentWindowId](KmerInfo& ki) { return ki.pos < currentWindowId; }
                    ),
                    heapWindow.end());
                std::make_heap(heapWindow.begin(), heapWindow.end(), KIHeap_cmp);
              }

              //Hash kmers
              hash_t hashFwd = CommonFunc::getHash(seq + i, kmerSize); 
              hash_t hashBwd;

              if(alphabetSize == 4)
                hashBwd = CommonFunc::getHash(seqRev.get() + (len - i - kmerSize), kmerSize);
              else  //proteins
                hashBwd = std::numeric_limits<hash_t>::max();   //Pick a dummy high value so that it is ignored later

              //Take minimum value of kmer and its reverse complement
              hash_t currentKmer = std::min(hashFwd, hashBwd);
              

              //Check the strand of this minimizer hash value
              auto currentStrand = hashFwd < hashBwd ? strnd::FWD : strnd::REV;

              //If front minimum is not in the current window, remove it
              if (!Q.empty() && std::get<2>(Q.front()) <  currentWindowId)
              {
                const auto [leaving_hash, leaving_strand, _] = Q.front();

                if (!sortedHashes.empty() && leaving_hash <= sortedHashes.back())
                {

                  auto& leaving_entry = windowMap.find(leaving_hash)->second;

                  // Check if this is the only occurence of this hash in the window
                  if (leaving_entry.count() == 1)
                  {
                    leaving_entry.mi.wpos_end = currentWindowId;
                    minmerIndex.push_back(leaving_entry.mi);
                    windowMap.erase(leaving_hash);
                    sortedHashes.erase(
                        std::lower_bound(sortedHashes.begin(), sortedHashes.end(), leaving_hash));
                  }
                  else
                  {
                    // Not removing hash, but need to adjust the strand
                    if (leaving_entry.mi.strand - leaving_strand == 0
                            || leaving_entry.mi.strand == 0)
                    {
                      leaving_entry.mi.wpos_end = currentWindowId;
                      minmerIndex.push_back(leaving_entry.mi);
                      leaving_entry.mi.wpos = currentWindowId;
                      leaving_entry.mi.wpos_end = -1;
                    }
                    leaving_entry.mi.strand -= leaving_strand;

                    // Remove position from poslist (FIFO head advance)
                    leaving_entry.occHead++;
                    if (leaving_entry.occHead == leaving_entry.occ.size()) {
                      leaving_entry.occ.clear();
                      leaving_entry.occHead = 0;
                    } else if (leaving_entry.occHead >= 32
                               && leaving_entry.occHead * 2 >= leaving_entry.occ.size()) {
                      leaving_entry.occ.erase(leaving_entry.occ.begin(),
                                              leaving_entry.occ.begin() + leaving_entry.occHead);
                      leaving_entry.occHead = 0;
                    }
                  }
                }
                Q.pop_front();
              }

              if (seq[i+kmerSize-1] == 'N')
              {
                ambig_kmer_count = kmerSize;
              }
              //Consider non-symmetric kmers only
              if(hashBwd != hashFwd && ambig_kmer_count == 0)
              {
                // Add current hash to window
                Q.push_back(std::make_tuple(currentKmer, currentStrand, i)); 

                // Check if current kmer is already in the map. Membership implies
                // currentKmer <= the max in-window hash, so skip the probe when
                // it is larger (the common case: the window keeps the smallest).
                auto kmer_it = (!sortedHashes.empty() && currentKmer <= sortedHashes.back())
                    ? windowMap.find(currentKmer) : windowMap.end();
                if (kmer_it != windowMap.end())
                {
                  auto& current_entry = kmer_it->second;
                  current_entry.occ.emplace_back(KmerInfo {currentKmer, seqCounter, i, currentStrand});
                  // Not removing hash, but need to adjust the strand
                  if (current_entry.mi.strand + currentStrand == 0
                          || current_entry.mi.strand == 0)
                  {
                    current_entry.mi.wpos_end = currentWindowId;
                    minmerIndex.push_back(current_entry.mi);
                    current_entry.mi.wpos = currentWindowId;
                    current_entry.mi.wpos_end = -1;
                  }
                  current_entry.mi.strand += currentStrand;
                }
                // Going in the heap
                else 
                {
                  heapWindow.emplace_back(KmerInfo {currentKmer, seqCounter, i, currentStrand});
                  std::push_heap(heapWindow.begin(), heapWindow.end(), KIHeap_cmp);
                }
              }
              if (ambig_kmer_count > 0)
              {
                ambig_kmer_count--;
              }
              



              // Add kmers from heap to window until full
              if(currentWindowId >= 0)
              {
                // Ignore expired kmers
                while (!heapWindow.empty() && heapWindow.front().pos < currentWindowId)
                {
                  std::pop_heap(heapWindow.begin(), heapWindow.end(), KIHeap_cmp);
                  heapWindow.pop_back(); 
                }

                //TODO leq?
                if (!sortedHashes.empty() && heapWindow.size() > 0
                    && sortedHashes.size() == (size_t)sketchSize
                    && (heapWindow.front().hash < sortedHashes.back()))
                {
                  auto largest_it = windowMap.find(sortedHashes.back());
                  auto& largest = largest_it->second;
                  // Add largest to index
                  largest.mi.wpos_end = currentWindowId;
                  minmerIndex.push_back(largest.mi);

                  // Add kmers back to heap
                  for (size_t oi = largest.occHead; oi < largest.occ.size(); ++oi)
                  {
                    if (largest.occ[oi].pos > currentWindowId) {
                        heapWindow.push_back(largest.occ[oi]);
                        std::push_heap(heapWindow.begin(), heapWindow.end(), KIHeap_cmp);
                    }
                  }

                  // Remove from window
                  windowMap.erase(largest_it);
                  sortedHashes.pop_back();
                }

                while (!heapWindow.empty() && sortedHashes.size() < (size_t)sketchSize)
                {
                  if (heapWindow.front().pos < currentWindowId)
                  {
                    std::pop_heap(heapWindow.begin(), heapWindow.end(), KIHeap_cmp);
                    heapWindow.pop_back();
                  }
                  // Add kmers of same value
                  const KmerInfo newKmer = heapWindow.front();
                  auto w_ins = windowMap.try_emplace(newKmer.hash);
                  auto& windowEntry = w_ins.first->second;
                  if (w_ins.second) {
                    sortedHashes.insert(
                        std::lower_bound(sortedHashes.begin(), sortedHashes.end(), newKmer.hash),
                        newKmer.hash);
                  }
                  windowEntry.mi = MinmerInfo{newKmer.hash, currentWindowId, -1, seqCounter, 0};
                  while (!heapWindow.empty() && heapWindow.front().hash == newKmer.hash)
                  {
                    windowEntry.occ.push_back(heapWindow.front());
                    windowEntry.mi.strand += heapWindow.front().strand;
                    std::pop_heap(heapWindow.begin(), heapWindow.end(), KIHeap_cmp);
                    heapWindow.pop_back();
                  }
                }
              }
            }

            // Add remaining open minmer windows (ascending hash order)
            uint64_t rank = 1;
            for (size_t si = 0; si < sortedHashes.size() && rank <= (uint64_t)sketchSize; ++si, ++rank)
            {
              auto& e = windowMap.find(sortedHashes[si])->second;
              if (e.mi.wpos != -1)
              {
                e.mi.wpos_end = len - kmerSize + 1;
                minmerIndex.push_back(e.mi);
              }
            }

            //// TODO Not sure why these are occuring but they are a bug
            minmerIndex.erase(
                std::remove_if(
                  minmerIndex.begin(), 
                  minmerIndex.end(), 
                  [](auto& mi) { return mi.wpos < 0 || mi.wpos_end < 0 || mi.wpos == mi.wpos_end; }),
                minmerIndex.end());


            //// Split up windows longer than windowSize into chunks of windowSize or less
            std::vector<MinmerInfo> chunkedMIs;
            std::for_each(minmerIndex.begin(), minmerIndex.end(), [&chunkedMIs, windowSize, kmerSize] (auto& mi) {
              mi.strand = mi.strand < 0 ? (mi.strand == 0 ? strnd::AMBIG : strnd::REV) : strnd::FWD;
              if (mi.wpos_end > mi.wpos + windowSize) {
                for (int chunk = 0; chunk < std::ceil(float(mi.wpos_end - mi.wpos) / float(windowSize)); chunk++) {
                  chunkedMIs.push_back(
                    MinmerInfo{
                      mi.hash, 
                      mi.wpos + chunk*windowSize, 
                      std::min(mi.wpos + chunk*windowSize + windowSize, mi.wpos_end),
                      mi.seqId, 
                      mi.strand
                    } 
                  );
                }
              }
            });
            minmerIndex.erase(
                std::remove_if(
                  minmerIndex.begin(), 
                  minmerIndex.end(), 
                  [windowSize](auto& mi) { return mi.wpos_end - mi.wpos > windowSize; }),
                minmerIndex.end());
            minmerIndex.insert(minmerIndex.end(), chunkedMIs.begin(), chunkedMIs.end());

            // Sort the index based on start position
            std::sort(minmerIndex.begin(), minmerIndex.end(), [](auto& l, auto& r) {return std::tie(l.wpos, l.wpos_end) < std::tie(r.wpos, r.wpos_end);});

            //// No duplicate windows
            //// TODO These should not be occurring. They happen rarely, so just deleting them for now
            //// but need to fix eventually 
            minmerIndex.erase(
                std::unique(
                  minmerIndex.begin(), 
                  minmerIndex.end(), 
                  [](auto& l, auto& r) { return (l.wpos == r.wpos) && (l.hash == r.hash); }),
                minmerIndex.end());

          }

        /**
          * @brief           Functor for comparing tuples by single index layer
          * @tparam layer    Tuple's index which is used for comparison
          * @tparam op       comparator, default as std::less
          */
        template<size_t layer, template<typename> class op = std::less>
        struct TpleComp {
            //Compare two tuples using their values
            template<typename T>
            bool operator()(T const &t1, T const &t2) {
                return op<typename std::tuple_element<layer, T>::type>()(std::get<layer>(t1), std::get<layer>(t2));
            }
        };

        /**
         * @brief                   computes the total size of reference in bytes
         * @param[in] refSequences  vector of reference files
         * @return                  total size
         */
        inline uint64_t getReferenceSize(const std::vector<std::string> &refSequences) {
            uint64_t count = 0;

            for (auto &f : refSequences) {
                //Open the file as binary, and set the position to end
                std::ifstream in(f, std::ifstream::ate | std::ifstream::binary);

                //the position of the current character
                count += (uint64_t) (in.tellg());
            }

            return count;
        }

        // Splitting
        template<typename Out>
        void split(const std::string &s, char delim, Out result) {
            std::stringstream ss(s);
            std::string item;
            while (std::getline(ss, item, delim)) {
                *(result++) = item;
            }
        }

        std::vector<std::string> split(const std::string &s, char delim) {
            std::vector<std::string> elems;
            split(s, delim, std::back_inserter(elems));
            return elems;
        }
    }
}

#endif
