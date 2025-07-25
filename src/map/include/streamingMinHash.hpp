/**
 * @file    streamingMinHash.hpp
 * @brief   Efficient streaming MinHash implementation with binary heap
 * @details Supports parallel sequence processing and group-wise sketch merging
 */

#ifndef STREAMING_MINHASH_HPP
#define STREAMING_MINHASH_HPP

#include <queue>
#include <vector>
#include <algorithm>
#include <mutex>
#include <unordered_map>
#include "map/include/base_types.hpp"
#include "map/include/map_parameters.hpp"
#include "map/include/sequenceIds.hpp"
#include "common/progress.hpp"

// Forward declarations
namespace skch {
    namespace CommonFunc {
        void makeUpperCaseAndValidDNA(char* seq, offset_t len);
        void reverseComplement(const char* src, char* dest, int length);
        hash_t getHash(const char* seq, int length);
    }
}

namespace skch {

/**
 * @class StreamingMinHash
 * @brief Maintains top-k minimum hash values using a max-heap
 */
class StreamingMinHash {
private:
    std::priority_queue<hash_t> maxHeap;  // Max-heap to keep k smallest values
    size_t sketchSize;
    size_t windowSize;
    std::mutex heapMutex;  // For thread-safe operations if needed

public:
    StreamingMinHash() : sketchSize(0), windowSize(0) {}  // Default constructor
    StreamingMinHash(size_t k, size_t w) : sketchSize(k), windowSize(w) {}
    
    // Copy constructor
    StreamingMinHash(const StreamingMinHash& other) 
        : maxHeap(other.maxHeap), sketchSize(other.sketchSize), windowSize(other.windowSize) {}
    
    // Move constructor
    StreamingMinHash(StreamingMinHash&& other) noexcept
        : maxHeap(std::move(other.maxHeap)), sketchSize(other.sketchSize), windowSize(other.windowSize) {}
    
    // Copy assignment
    StreamingMinHash& operator=(const StreamingMinHash& other) {
        if (this != &other) {
            maxHeap = other.maxHeap;
            sketchSize = other.sketchSize;
            windowSize = other.windowSize;
        }
        return *this;
    }
    
    // Move assignment
    StreamingMinHash& operator=(StreamingMinHash&& other) noexcept {
        if (this != &other) {
            maxHeap = std::move(other.maxHeap);
            sketchSize = other.sketchSize;
            windowSize = other.windowSize;
        }
        return *this;
    }

    /**
     * @brief Add a hash value to the sketch (thread-safe)
     * @param hash The hash value to add
     * @return true if the hash was added (is in top-k), false otherwise
     */
    bool add(hash_t hash) {
        std::lock_guard<std::mutex> lock(heapMutex);
        return add_unsafe(hash);
    }
    
    /**
     * @brief Add a hash value to the sketch (NOT thread-safe, for single-threaded use)
     * @param hash The hash value to add
     * @return true if the hash was added (is in top-k), false otherwise
     */
    bool add_unsafe(hash_t hash) {
        if (maxHeap.size() < sketchSize) {
            maxHeap.push(hash);
            return true;
        } else if (hash < maxHeap.top()) {
            maxHeap.pop();
            maxHeap.push(hash);
            return true;
        }
        return false;
    }

    /**
     * @brief Get the current sketch as a sorted vector
     * @return Vector of hash values in ascending order
     */
    std::vector<hash_t> getSketch() const {
        std::vector<hash_t> result;
        result.reserve(maxHeap.size());
        
        // Copy heap to preserve it
        auto heapCopy = maxHeap;
        while (!heapCopy.empty()) {
            result.push_back(heapCopy.top());
            heapCopy.pop();
        }
        
        // Sort in ascending order
        std::sort(result.begin(), result.end());
        return result;
    }

    /**
     * @brief Merge another sketch into this one
     * @param other The sketch to merge
     */
    void merge(const StreamingMinHash& other) {
        auto otherSketch = other.getSketch();
        for (hash_t hash : otherSketch) {
            add(hash);
        }
    }

    size_t size() const { return maxHeap.size(); }
    bool empty() const { return maxHeap.empty(); }
    void clear() { 
        std::priority_queue<hash_t> empty;
        std::swap(maxHeap, empty);
    }
};

/**
 * @class GroupedStreamingMinHash
 * @brief Manages MinHash sketches grouped by sequence prefix
 */
class GroupedStreamingMinHash {
private:
    std::unordered_map<int, StreamingMinHash> groupSketches;
    size_t sketchSize;
    size_t windowSize;
    std::mutex groupMutex;

public:
    GroupedStreamingMinHash(size_t k, size_t w) : sketchSize(k), windowSize(w) {}

    /**
     * @brief Process a sequence and update the appropriate group sketch
     * @param seq The sequence data
     * @param len Sequence length
     * @param seqId Sequence ID
     * @param groupId Group ID (based on prefix)
     * @param kmerSize K-mer size
     * @param alphabetSize Alphabet size (4 for DNA)
     * @param progress Progress meter for updates
     */
    void processSequence(const char* seq, offset_t len, seqno_t seqId, 
                        int groupId, int kmerSize, int alphabetSize,
                        progress_meter::ProgressMeter* progress = nullptr) {
        
        // Ensure group exists
        {
            std::lock_guard<std::mutex> lock(groupMutex);
            if (groupSketches.find(groupId) == groupSketches.end()) {
                groupSketches.emplace(groupId, StreamingMinHash(sketchSize, windowSize));
            }
        }

        // Buffer for current k-mer (small, reusable)
        std::vector<char> kmerBuf(kmerSize);
        std::vector<char> seqRev(kmerSize);

        // Track ambiguous k-mers
        int ambig_kmer_count = 0;
        
        // Check initial k-mer for N's
        for (int j = 0; j < kmerSize && j < len; j++) {
            char c = std::toupper(seq[j]);
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
                ambig_kmer_count = kmerSize;
                break;
            }
        }

        // Process all k-mers in the sequence
        for (offset_t i = 0; i <= len - kmerSize; i++) {
            if (progress) progress->increment(1);

            // Check for N at the end of current k-mer
            char endChar = std::toupper(seq[i + kmerSize - 1]);
            if (endChar != 'A' && endChar != 'C' && endChar != 'G' && endChar != 'T') {
                ambig_kmer_count = kmerSize;
            }

            if (ambig_kmer_count == 0) {
                // Copy and uppercase k-mer
                for (int j = 0; j < kmerSize; j++) {
                    kmerBuf[j] = std::toupper(seq[i + j]);
                }
                
                // Hash forward k-mer
                hash_t hashFwd = CommonFunc::getHash(kmerBuf.data(), kmerSize);
                hash_t hashBwd = std::numeric_limits<hash_t>::max();

                // Hash reverse complement for DNA
                if (alphabetSize == 4) {
                    CommonFunc::reverseComplement(kmerBuf.data(), seqRev.data(), kmerSize);
                    hashBwd = CommonFunc::getHash(seqRev.data(), kmerSize);
                }

                // Use canonical k-mer (minimum of forward and reverse)
                if (hashFwd != hashBwd) {
                    hash_t canonicalHash = std::min(hashFwd, hashBwd);
                    
                    // Add to group sketch
                    std::lock_guard<std::mutex> lock(groupMutex);
                    groupSketches[groupId].add(canonicalHash);
                }
            }

            if (ambig_kmer_count > 0) {
                ambig_kmer_count--;
            }
        }
    }

    /**
     * @brief Get MinmerInfo objects for a specific group
     * @param groupId The group ID
     * @param seqId Sequence ID to use for MinmerInfo
     * @return Vector of MinmerInfo objects
     */
    std::vector<MinmerInfo> getGroupMinmers(int groupId, seqno_t seqId) const {
        std::vector<MinmerInfo> result;
        
        auto it = groupSketches.find(groupId);
        if (it != groupSketches.end()) {
            auto sketch = it->second.getSketch();
            
            // Convert hash values to MinmerInfo objects
            // Note: We're creating synthetic window positions since we don't track them
            offset_t pos = 0;
            for (hash_t hash : sketch) {
                result.push_back(MinmerInfo{
                    hash,
                    pos,              // wpos
                    static_cast<offset_t>(pos + windowSize), // wpos_end
                    seqId,
                    strnd::FWD       // Default to forward strand
                });
                pos += windowSize;
            }
        }
        
        return result;
    }

    /**
     * @brief Get all group sketches
     * @return Map of group ID to sketch
     */
    std::unordered_map<int, std::vector<hash_t>> getAllGroupSketches() const {
        std::unordered_map<int, std::vector<hash_t>> result;
        
        for (const auto& [groupId, sketch] : groupSketches) {
            result[groupId] = sketch.getSketch();
        }
        
        return result;
    }

    /**
     * @brief Clear all sketches
     */
    void clear() {
        std::lock_guard<std::mutex> lock(groupMutex);
        groupSketches.clear();
    }
};

/**
 * @brief Parallel sequence processor using streaming MinHash
 * @param sequences Vector of sequence data with metadata
 * @param param Sketch parameters
 * @param idManager Sequence ID manager for group assignment
 * @param progress Progress meter
 * @return Populated minmer index
 */
template<typename SeqContainer>
std::vector<MinmerInfo> parallelStreamingMinHash(
    const std::vector<SeqContainer>& sequences,
    const Parameters& param,
    const skch::SequenceIdManager& idManager,
    progress_meter::ProgressMeter* progress = nullptr) {
    
    GroupedStreamingMinHash groupedSketch(param.sketchSize, param.windowLength);
    
    // Process sequences in parallel
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < sequences.size(); i++) {
        const auto& seqContainer = sequences[i];
        int groupId = idManager.getRefGroup(seqContainer.seqId);
        
        groupedSketch.processSequence(
            seqContainer.seq.data(),
            seqContainer.seq.length(),
            seqContainer.seqId,
            groupId,
            param.kmerSize,
            param.alphabetSize,
            progress
        );
    }
    
    // Collect all minmers from all groups
    std::vector<MinmerInfo> allMinmers;
    for (const auto& seqContainer : sequences) {
        int groupId = idManager.getRefGroup(seqContainer.seqId);
        auto groupMinmers = groupedSketch.getGroupMinmers(groupId, seqContainer.seqId);
        allMinmers.insert(allMinmers.end(), groupMinmers.begin(), groupMinmers.end());
    }
    
    return allMinmers;
}

} // namespace skch

#endif // STREAMING_MINHASH_HPP