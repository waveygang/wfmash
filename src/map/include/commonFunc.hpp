/**
 * @file    commonFunc.hpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef COMMON_FUNC_HPP
#define COMMON_FUNC_HPP

#include <vector>
#include <algorithm>
#include <deque>
#include <cmath>
#include <fstream>
#include <limits>

//Own includes
#include "map/include/map_parameters.hpp"

//External includes
#include "common/murmur3.h"
#include "common/prettyprint.hpp"

#include "common/wflign/src/rkmh.hpp"

namespace skch {
    /**
     * @namespace skch::CommonFunc
     * @brief     Implements frequently used common functions
     */
    namespace CommonFunc {
        //seed for murmerhash
        const int seed = 42;

        /**
         * @brief   reverse complement of kmer (borrowed from mash)
         * @note    assumes dest is pre-allocated
         */
        inline void reverseComplement(const char *src, char *dest, int length) {
            for (int i = 0; i < length; i++) {
                char base = src[i];

                switch (base) {
                    case 'A':
                        base = 'T';
                        break;
                    case 'C':
                        base = 'G';
                        break;
                    case 'G':
                        base = 'C';
                        break;
                    case 'T':
                        base = 'A';
                        break;
                    default:
                        break;
                }

                dest[length - i - 1] = base;
            }
        }

        /**
     * @brief               convert DNA or AA alphabets to upper case, converting non-canonical DNA bases to N
     * @param[in]   seq     pointer to input sequence
     * @param[in]   len     length of input sequence
     */
        inline void makeUpperCaseAndValidDNA(char *seq, offset_t len) {
            for (int i = 0; i < len; i++) {
                if (seq[i] > 96 && seq[i] < 123) {
                    seq[i] -= 32;
                }

                if (rkmh::valid_dna[seq[i]]) {
                    seq[i] = 'N';
                }
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
         * @brief       compute winnowed minimizers from a given sequence and add to the index
         * @param[out]  minimizerIndex  minimizer table storing minimizers and their position as we compute them
         * @param[in]   seq             pointer to input sequence
         * @param[in]   len             length of input sequence
         * @param[in]   kmerSize
         * @param[in]   windowSize
         * @param[in]   seqCounter      current sequence number, used while saving the position of minimizer
         */
        template<typename T>
        inline void addMinimizers(std::vector<T> &minimizerIndex,
                                  char *seq, offset_t len,
                                  int kmerSize,
                                  int windowSize,
                                  int alphabetSize,
                                  seqno_t seqCounter
                                  //const std::unordered_set<std::string>& high_freq_kmers
                                  ) {
            /**
             * Double-ended queue (saves minimum at front end)
             * Saves pair of the minimizer and the position of hashed kmer in the sequence
             * Position of kmer is required to discard kmers that fall out of current window
             */
            std::deque<std::pair<MinimizerInfo, offset_t> > Q;

            makeUpperCaseAndValidDNA(seq, len);

            //Compute reverse complement of seq
            char *seqRev = new char[len];

            if (alphabetSize == 4) //not protein
                CommonFunc::reverseComplement(seq, seqRev, len);

            for (offset_t i = 0; i < len - kmerSize + 1; i++) {
                //The serial number of current sliding window
                //First valid window appears when i = windowSize - 1
                offset_t currentWindowId = i - windowSize + 1;

                //Hash kmers
                hash_t hashFwd = CommonFunc::getHash(seq + i, kmerSize);
                hash_t hashBwd;

                if (alphabetSize == 4)
                    hashBwd = CommonFunc::getHash(seqRev + len - i - kmerSize, kmerSize);
                else  //proteins
                    hashBwd = std::numeric_limits<hash_t>::max();   //Pick a dummy high value so that it is ignored later

//#define DEBUG_WINNOWING
#ifdef DEBUG_WINNOWING
                std::cout << "pos: " << i << std::endl;
                std::cout << "kmers: ";
                for (uint64_t j = 0; j < kmerSize; ++j) {
                    std::cout << seq[i + j];
                }
                std::cout << " --> " << hashFwd << " - " << hashBwd << std::endl;

                std::cout << "Q1" << std::endl;
                for(auto iter = Q.begin(); iter != Q.end(); ++iter) {
                    std::cout << iter->second << " " << " " << iter->first.hash <<  " " << iter->first.wpos << std::endl;
                }
                std::cout << std::endl;
#endif

                //Consider non-symmetric kmers only
                if (hashBwd != hashFwd) {
                    //Take minimum value of kmer and its reverse complement
                    hash_t currentKmer = std::min(hashFwd, hashBwd);

                    /*double order = (hashFwd < hashBwd) ?
                                   applyWeight(seq + i, kmerSize, hashFwd, high_freq_kmers) :
                                   applyWeight(seqRev + len - i - kmerSize, kmerSize, hashBwd, high_freq_kmers);*/

                    //Hashes less than equal to currentKmer are not required
                    //Remove them from Q (back)
                    //while (!Q.empty() && Q.back().first.order > order)
                    while (!Q.empty() && Q.back().first.hash > currentKmer)
                        Q.pop_back();

#ifdef DEBUG_WINNOWING
                    std::cout << "Q2" << std::endl;
                    for(auto iter = Q.begin(); iter != Q.end(); ++iter) {
                        std::cout << iter->second << " " << " " << iter->first.hash <<  " " << iter->first.wpos << std::endl;
                    }
                    std::cout << std::endl;
#endif

                    //Check the strand of this minimizer hash value
                    auto currentStrand = hashFwd < hashBwd ? strnd::FWD : strnd::REV;

                    //Push currentKmer and position to back of the queue
                    //-1 indicates the dummy window # (will be updated later)
                    Q.push_back(std::make_pair(
                            //MinimizerInfo{currentKmer, seqCounter, -1, currentStrand, order},
                            MinimizerInfo{currentKmer, seqCounter, -1, currentStrand},
                            i));

#ifdef DEBUG_WINNOWING
                    std::cout << "Q3" << std::endl;
                    for(auto iter = Q.begin(); iter != Q.end(); ++iter) {
                        std::cout << iter->second << " " << " " << iter->first.hash <<  " " << iter->first.wpos << std::endl;
                    }
                    std::cout << std::endl;
#endif

                    //If front minimum is not in the current window, remove it
                    if (!Q.empty() && Q.front().second <= i - windowSize) {
                        while (!Q.empty() && Q.front().second <= i - windowSize)
                            Q.pop_front();
#ifdef DEBUG_WINNOWING
                        std::cout << "Q4" << std::endl;
                        for(auto iter = Q.begin(); iter != Q.end(); ++iter) {
                            std::cout << iter->second << " " << " " << iter->first.hash <<  " " << iter->first.wpos << std::endl;
                        }
                        std::cout << std::endl;
#endif

                        // Robust-winnowing
                        //while (Q.size() > 1 && Q.begin()->first.order == (++Q.begin())->first.order)
                        while (Q.size() > 1 && Q.begin()->first.hash == (++Q.begin())->first.hash)
                            Q.pop_front();
                    }

#ifdef DEBUG_WINNOWING
                    std::cout << "Q5" << std::endl;
                    for(auto iter = Q.begin(); iter != Q.end(); ++iter) {
                        std::cout << iter->second << " " << " " << iter->first.hash <<  " " << iter->first.wpos << std::endl;
                    }
                    std::cout << std::endl;
#endif

                    //Select the minimizer from Q and put into index
                    if (currentWindowId >= 0) {
                        //We save the minimizer if we are seeing it for first time
                        if (minimizerIndex.empty() || minimizerIndex.back() != Q.front().first) {
                            //Update the window position in this minimizer
                            //This step also ensures we don't re-insert the same minimizer again
                            Q.front().first.wpos = currentWindowId;
                            minimizerIndex.push_back(Q.front().first);

#ifdef DEBUG_WINNOWING
                            std::cout << "PUSHED: " << Q.front().first.wpos << " " << Q.front().first.hash << std::endl;
#endif
                        }
                    }

#ifdef DEBUG_WINNOWING
                    std::cout << "Q - FINAL" << std::endl;
                    for(auto iter = Q.begin(); iter != Q.end(); ++iter) {
                        std::cout << iter->second << " " << " " << iter->first.hash <<  " " << iter->first.wpos << std::endl;
                    }
                    std::cout << std::endl;

                    std::cout << "minimizerIndex" << std::endl;
                    for(auto iter = minimizerIndex.begin(); iter != minimizerIndex.end(); ++iter) {
                        std::cout << iter->wpos << " " << iter->hash << std::endl;
                    }
                    std::cout << std::endl;
#endif
                }
#ifdef DEBUG_WINNOWING
                std::cout << "--------------------------------------------------------" << std::endl;
#endif
            }

#ifdef DEBUG
            std::cerr << "INFO, skch::CommonFunc::addMinimizers, inserted minimizers for sequence id = " << seqCounter << "\n";
#endif

            delete[] seqRev;
        }

        /**
         * @brief       compute winnowed minimizers from a given sequence and add to the index
         * @param[out]  minimizerIndex  minimizer table storing minimizers and their position as we compute them
         * @param[in]   seq             pointer to input sequence
         * @param[in]   len             length of input sequence
         * @param[in]   kmerSize
         * @param[in]   samplingFactor
         * @param[in]   seqCounter      current sequence number, used while saving the position of minimizer
         */
        template<typename T>
        inline void addWorldMinimizers(std::vector<T> &minimizerIndex,
                                       char *seq, offset_t len,
                                       int kmerSize,
                                       int samplingFactor,
                                       int alphabetSize,
                                       seqno_t seqCounter) {

            makeUpperCaseAndValidDNA(seq, len);

            //Compute reverse complement of seq
            char *seqRev = new char[len];

            // get our sampling fraction
            hash_t samplingBound = std::numeric_limits<hash_t>::max() / samplingFactor;

            if (alphabetSize == 4) //not protein
                CommonFunc::reverseComplement(seq, seqRev, len);

            for (offset_t i = 0; i < len - kmerSize + 1; i++) {
                //Hash kmers
                hash_t hashFwd = CommonFunc::getHash(seq + i, kmerSize);
                hash_t hashBwd;

                if (alphabetSize == 4)
                    hashBwd = CommonFunc::getHash(seqRev + len - i - kmerSize, kmerSize);
                else  //proteins
                    hashBwd = std::numeric_limits<hash_t>::max();   //Pick a dummy high value so that it is ignored later
                if (hashBwd != hashFwd) { // consider non-symmetric kmers only
                    if (hashFwd < samplingBound) {
                        minimizerIndex.push_back(MinimizerInfo{hashFwd, seqCounter, i, strnd::FWD});
                    }
                    if (hashBwd < samplingBound) {
                        minimizerIndex.push_back(MinimizerInfo{hashBwd, seqCounter, i, strnd::REV});
                    }
                }
            }
            delete[] seqRev;
        }


        /**
         * @brief       compute winnowed minimizers from a given sequence and add to the index
         * @param[out]  minimizerIndex  minimizer table storing minimizers and their position as we compute them
         * @param[in]   seq             pointer to input sequence
         * @param[in]   len             length of input sequence
         * @param[in]   kmerSize
         * @param[in]   samplingFactor
         * @param[in]   seqCounter      current sequence number, used while saving the position of minimizer
         * @param[in]   spaced_seeds    A vector of spaced seeds from ALeS
         */
        template<typename T>
        inline void addSpacedSeedWorldMinimizers(std::vector<T> &minimizerIndex,
                                                 char *seq,
                                                 offset_t len,
                                                 int kmerSize,
                                                 int samplingFactor,
                                                 int alphabetSize,
                                                 seqno_t seqCounter,
                                                 const std::vector<ales::spaced_seed> &spaced_seeds) {

            makeUpperCaseAndValidDNA(seq, len);
            size_t minimizer_range_start = minimizerIndex.size();
            
            //Compute reverse complement of seq
            char* seqRev = new char[len];

            // get our sampling fraction
            hash_t samplingBound = std::numeric_limits<hash_t>::max() / samplingFactor;

            // not protein
            if (alphabetSize == 4) {
                CommonFunc::reverseComplement(seq, seqRev, len);
            }

            // we increment this value once a spaced seed would go out of bounds of the sequence
            int handled_seed_count = 0;
            // start position on the sequence
            offset_t i = 0;
            // while we still have spaced seeds whose hashes can be computers
            while (handled_seed_count < spaced_seeds.size()) {

                for (uint32_t spaced_seed_number = 0; spaced_seed_number < spaced_seeds.size(); spaced_seed_number++) {
                    const auto &s = spaced_seeds[spaced_seed_number];
                    size_t seed_length = s.length;
                    char *ss = s.seed;

                    // if the (spaced) seed would go out of bounds, skip it
                    if (i + seed_length >= len) {
                        handled_seed_count = handled_seed_count + 1;
                        continue;
                    }

                    char *forward_start_char = seq + i;
                    char *reverse_start_char = seqRev + len - i - seed_length;
                    char new_forward_kmer[seed_length];
                    char new_reverse_kmer[seed_length];

                    for (size_t j = 0; j < seed_length; ++j, ++ss, ++forward_start_char) {
                        new_forward_kmer[j] = *ss == '1' ? *forward_start_char : '*';
                    }
                    ss = s.seed;
                    for (size_t j = 0; j < seed_length; ++j, ++ss, ++reverse_start_char) {
                        new_reverse_kmer[j] = *ss == '1' ? *reverse_start_char : '*';
                    }
                    ss = s.seed;// reset the seed for the next iteration of the loop

                    /* debug print
                       std::cerr << seed_length << " " << s.seed
                       << " forward " << new_forward_kmer
                       << " reverse " << new_reverse_kmer << std::endl;
                    */


                    //Hash kmers
                    hash_t hashFwd = CommonFunc::getHash(&new_forward_kmer[0], seed_length);
                    hash_t hashBwd;

                    if (alphabetSize == 4)
                        hashBwd = CommonFunc::getHash(&new_reverse_kmer[0], seed_length);
                    else                                             //proteins
                        hashBwd = std::numeric_limits<hash_t>::max();//Pick a dummy high value so that it is ignored later

                    // Consider non-symmetric kmers only
                    if (hashBwd != hashFwd) {
                        if (hashFwd < samplingBound) {
                            minimizerIndex.push_back(MinimizerInfo{hashFwd, seqCounter, i, strnd::FWD});
                        }
                        if (hashBwd < samplingBound) {
                            minimizerIndex.push_back(MinimizerInfo{hashBwd, seqCounter, i, strnd::REV});
                        }
                    }
                }

                i = i + 1;
            }

            delete[] seqRev;
        }

        /**
         * @brief       compute winnowed minimizers from a given sequence and add to the index using spaced seeds
         * @param[out]  minimizerIndex  minimizer table storing minimizers and their position as we compute them
         * @param[in]   seq             pointer to input sequence
         * @param[in]   len             length of input sequence
         * @param[in]   kmerSize
         * @param[in]   windowSize
         * @param[in]   seqCounter      current sequence number, used while saving the position of minimizer
         * @param[in]   spaced_seeds    A vector of spaced seeds from ALeS
         */
        template <typename T>
        void addSpacedSeedMinimizers(std::vector<T> &minimizerIndex,
                                     char* seq,
                                     offset_t len,
                                     int kmerSize,
                                     int windowSize,
                                     int alphabetSize,
                                     seqno_t seqCounter,
                                     const std::vector<ales::spaced_seed>& spaced_seeds
                                     )
        {

          makeUpperCaseAndValidDNA(seq, len);
          size_t minimizer_range_start = minimizerIndex.size();

          //Compute reverse complement of seq
          char* seqRev = new char[len];

          auto extract_kmer = [](char* thing, size_t len) {
            std::string the_string;
            for (size_t i=0; i<len; i++, thing++)
              the_string.push_back(*thing);

            return the_string;
          };

          if(alphabetSize == 4) //not protein
            CommonFunc::reverseComplement(seq, seqRev, len);

          for (uint32_t spaced_seed_number=0; spaced_seed_number < spaced_seeds.size(); spaced_seed_number++) {
            /**
             * Double-ended queue (saves minimum at front end)
             * Saves pair of the minimizer and the position of hashed kmer in the sequence
             * Position of kmer is required to discard kmers that fall out of current window
             */
            std::deque< std::pair<MinimizerInfo, offset_t> > Q;
            const auto& s = spaced_seeds[spaced_seed_number];
            size_t seed_length =  s.length;
            char* ss = s.seed;

            for (offset_t i = 0; i < len - seed_length + 1; i++) {
              char* forward_start_char = seq+i;
              char* reverse_start_char = seqRev + len - i - seed_length;
              char new_forward_kmer[seed_length];
              char new_reverse_kmer[seed_length];

              for (size_t j=0; j<seed_length; ++j, ++ss, ++forward_start_char) {
                  new_forward_kmer[j] = *ss == '1' ? *forward_start_char : '*';
              }
              ss = s.seed;
              for (size_t j=0; j<seed_length; ++j, ++ss, ++reverse_start_char) {
                  new_reverse_kmer[j] = *ss == '1' ? *reverse_start_char : '*';
              }
              ss = s.seed; // reset the seed for the next iteration of the loop

              /* debug print
                 std::cerr << seed_length << " " << s.seed
                 << " forward " << extract_kmer(seq+i, seed_length) << " " << new_forward_kmer
                 << " reverse " << extract_kmer(seqRev + len - i - seed_length, seed_length) << " " << new_reverse_kmer << std::endl;
              */

              offset_t currentWindowId = i - windowSize + 1;

              //Hash kmers
              hash_t hashFwd = CommonFunc::getHash(&new_forward_kmer[0], seed_length);
              hash_t hashBwd;

              if(alphabetSize == 4)
                hashBwd = CommonFunc::getHash(&new_reverse_kmer[0], seed_length);
              else  //proteins
                hashBwd = std::numeric_limits<hash_t>::max();   //Pick a dummy high value so that it is ignored later

              // Consider non-symmetric kmers only
              if(hashBwd != hashFwd) {
                //Take minimum value of kmer and its reverse complement
                hash_t currentKmer = std::min(hashFwd, hashBwd);

                //Check the strand of this minimizer hash value
                auto currentStrand = hashFwd < hashBwd ? strnd::FWD : strnd::REV;

                //If front minimum is not in the current window, remove it
                while(!Q.empty() && Q.front().second <=  i - windowSize)
                  Q.pop_front();

                // Hashes less than equal to currentKmer are not required
                // Remove them from Q (back)
                while(!Q.empty() && Q.back().first.hash >= currentKmer)
                  Q.pop_back();

                // Push currentKmer and position to back of the queue
                // -1 indicates the dummy window # (will be updated later)
                Q.push_back(std::make_pair(MinimizerInfo {currentKmer,seqCounter, -1, currentStrand}, i));

                // Select the minimizer from Q and put into index
                if(currentWindowId >= 0) {
                    //We save the minimizer if we are seeing it for first time
                    if(minimizerIndex.empty() || minimizerIndex.back() != Q.front().first)
                      {
                        //Update the window position in this minimizer
                        //This step also ensures we don't re-insert the same minimizer again
                        Q.front().first.wpos = currentWindowId;
                        minimizerIndex.push_back(Q.front().first);
                      }
                  }
              }
            }
          }

          // sort our minimizerIndex by window position
          std::sort(minimizerIndex.begin() + minimizer_range_start,
                    minimizerIndex.end(),
                    [](const MinimizerInfo& a, const MinimizerInfo& b) {
                        return a.wpos < b.wpos;
                    });

#ifdef DEBUG
          std::cout << "INFO, skch::CommonFunc::addMinimizers, inserted minimizers for sequence id = " << seqCounter << "\n";
#endif

          delete [] seqRev;
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
