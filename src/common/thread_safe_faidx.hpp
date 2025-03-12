#ifndef THREAD_SAFE_FAIDX_HPP
#define THREAD_SAFE_FAIDX_HPP

// Standard C++ includes
#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <cstring>
#include <atomic>
#include <sstream>
#include <cstdio>
#include <iostream>

// HTSlib includes
#include <htslib/bgzf.h>
#include <htslib/faidx.h>

namespace ts_faidx {

/**
 * @brief Format of the FASTA/FASTQ file
 */
enum class FileFormat {
    UNKNOWN,
    FASTA,
    FASTQ
};

/**
 * @brief Entry for a single sequence in the FASTA/FASTQ index
 * 
 * This corresponds to one line in the .fai file
 */
class IndexEntry {
public:
    std::string name;        // Sequence name
    int64_t length;          // Sequence length
    int64_t offset;          // Offset in the file to the first base
    int32_t line_bases;      // Number of bases on each line (excluding newline)
    int32_t line_width;      // Total width of each line including newlines
    int64_t qual_offset;     // Offset to quality values (FASTQ only)
    
    /**
     * @brief Calculate file offset for a given sequence position
     * 
     * @param position 0-based position in the sequence
     * @return int64_t File offset for this position
     */
    int64_t calculate_offset(int64_t position) const {
        // Calculate the offset in the file for a specific sequence position
        int64_t line_number = position / line_bases;
        int64_t line_position = position % line_bases;
        return offset + line_number * line_width + line_position;
    }
};

/**
 * @brief BGZF file handler that provides random access to compressed files
 * 
 * This class uses the BGZF implementation from htslib which provides
 * efficient random access to bgzip-compressed files
 */
class BGZFReader {
private:
    BGZF* bgzf_ = nullptr;
    std::string filename_;
    bool is_compressed_ = false;
    
public:
    /**
     * @brief Construct a new BGZFReader object
     * 
     * @param filename Path to the file
     * @param mode Open mode ("r" for read, "w" for write)
     */
    BGZFReader(const std::string& filename, const char* mode = "r") : filename_(filename) {
        // Open the file with BGZF
        bgzf_ = bgzf_open(filename.c_str(), mode);
        if (!bgzf_) {
            throw std::runtime_error("Failed to open file: " + filename);
        }
        
        // Check if it's compressed
        is_compressed_ = bgzf_compression(bgzf_) > 0;
    }
    
    /**
     * @brief Check if file is compressed
     * 
     * @return true if compressed
     */
    bool is_compressed() const {
        return is_compressed_;
    }
    
    /**
     * @brief Seek to position in the file
     * 
     * @param position File offset
     * @return true if successful
     */
    bool seek(int64_t position) {
        return bgzf_seek(bgzf_, position, SEEK_SET) >= 0;
    }
    
    /**
     * @brief Read data from the file
     * 
     * @param buffer Buffer to read into
     * @param length Number of bytes to read
     * @return int64_t Number of bytes read
     */
    int64_t read(void* buffer, size_t length) {
        return bgzf_read(bgzf_, buffer, length);
    }
    
    /**
     * @brief Close the file
     */
    void close() {
        if (bgzf_) {
            bgzf_close(bgzf_);
            bgzf_ = nullptr;
        }
    }
    
    /**
     * @brief Get the BGZF handle
     * 
     * @return BGZF* Pointer to the BGZF handle
     */
    BGZF* get_handle() {
        return bgzf_;
    }
    
    /**
     * @brief Set the cache size
     * 
     * @param cache_size Cache size in bytes
     */
    void set_cache_size(int cache_size) {
        bgzf_set_cache_size(bgzf_, cache_size);
    }
    
    /**
     * @brief Destructor closes the file
     */
    ~BGZFReader() {
        close();
    }
};

/**
 * @brief Parse a region string into sequence name and coordinates
 * 
 * @param region Region string (e.g., "chr1:1000-2000")
 * @param seq_name Output parameter for sequence name
 * @param start Output parameter for start position (0-based)
 * @param end Output parameter for end position (0-based, inclusive)
 * @return true if successfully parsed
 */
bool parse_region(const std::string& region, 
                 std::string& seq_name, 
                 int64_t& start, 
                 int64_t& end) {
    // Default values
    start = 0;
    end = INT64_MAX;
    
    // Find the colon separating sequence name from coordinates
    size_t colon_pos = region.find(':');
    
    if (colon_pos == std::string::npos) {
        // No coordinates, just sequence name
        seq_name = region;
        return true;
    }
    
    // Extract sequence name
    seq_name = region.substr(0, colon_pos);
    
    // Find the dash separating start and end
    size_t dash_pos = region.find('-', colon_pos);
    
    if (dash_pos == std::string::npos) {
        // No end coordinate, just start
        try {
            start = std::stoll(region.substr(colon_pos + 1));
            end = start + 1; // Single base
        } catch (const std::exception&) {
            return false;
        }
    } else {
        // Both start and end coordinates
        try {
            start = std::stoll(region.substr(colon_pos + 1, dash_pos - colon_pos - 1));
            end = std::stoll(region.substr(dash_pos + 1));
        } catch (const std::exception&) {
            return false;
        }
    }
    
    // Convert to 0-based coordinates if they are 1-based
    if (start > 0) start--; // Assuming 1-based input, convert to 0-based
    
    return true;
}

/**
 * @brief Thread-safe FASTA/FASTQ index reader
 * 
 * This class provides thread-safe access to FASTA/FASTQ files with random access
 * to specific regions. Each thread automatically gets its own file handle.
 */
class FastaReader {
private:
    std::string filename_;                    // Path to the FASTA/FASTQ file
    std::unordered_map<std::string, IndexEntry> entries_; // Sequence entries
    std::vector<std::string> sequence_names_; // List of sequence names in order
    FileFormat format_ = FileFormat::UNKNOWN; // Format of the file
    mutable std::mutex resource_mutex_;      // Mutex for thread safety
    
    // Thread-local storage for file handles
    mutable std::unordered_map<std::thread::id, std::unique_ptr<BGZFReader>> thread_local_files_;
    
    /**
     * @brief Get a thread-local file handle
     * 
     * @return BGZFReader* Pointer to the thread's file handle
     */
    BGZFReader* get_thread_local_file() const {
        std::thread::id this_id = std::this_thread::get_id();
        
        // Lock to safely check/modify the map
        std::lock_guard<std::mutex> lock(resource_mutex_);
        
        auto it = thread_local_files_.find(this_id);
        if (it == thread_local_files_.end()) {
            // Create a new file handle for this thread
            std::unique_ptr<BGZFReader> file_handle(new BGZFReader(filename_));
            BGZFReader* result = file_handle.get();
            thread_local_files_[this_id] = std::move(file_handle);
            return result;
        }
        
        return it->second.get();
    }
    
    /**
     * @brief Retrieve sequence or quality from file
     * 
     * @param entry Index entry for the sequence
     * @param offset File offset for the data
     * @param start Start position (0-based)
     * @param end End position (0-based, exclusive)
     * @return std::string The requested data
     */
    std::string retrieve_data(const IndexEntry& entry, 
                             int64_t offset, 
                             int64_t start, 
                             int64_t end) const {
        // Use htslib's faidx directly for reliable sequence fetching
        std::string region = entry.name + ":" + 
                            std::to_string(start+1) + "-" + 
                            std::to_string(end);
        
        faidx_t* fai = fai_load(filename_.c_str());
        if (!fai) {
            throw std::runtime_error("Failed to load index for " + filename_);
        }
        
        hts_pos_t seq_len;
        char* seq = fai_fetch64(fai, region.c_str(), &seq_len);
        
        if (!seq) {
            fai_destroy(fai);
            throw std::runtime_error("Failed to fetch sequence for region: " + region);
        }
        
        std::string result(seq, seq_len);
        free(seq);
        fai_destroy(fai);
        
        return result;
    }
    
    /**
     * @brief Build an index for the FASTA/FASTQ file
     * 
     * @param out_path Output path for the index file
     * @return true if successful
     */
    bool build_index_internal(const std::string& out_path) {
        // Open the input file
        std::ifstream input(filename_, std::ios::binary);
        if (!input.is_open()) {
            return false;
        }
        
        // Open the output file
        std::ofstream output(out_path);
        if (!output.is_open()) {
            return false;
        }
        
        // Variables for index building
        std::string line;
        std::string name;
        int64_t offset = 0;
        int64_t seq_len = 0;
        int64_t seq_offset = 0;
        int64_t qual_offset = 0;
        int32_t line_len = 0;
        int32_t line_bases = 0;
        
        // State tracking
        enum State { NONE, IN_SEQ, IN_QUAL };
        State state = NONE;
        
        // Scan the file
        while (std::getline(input, line)) {
            int64_t line_bytes = line.length() + 1; // +1 for newline
            
            if (line.empty()) {
                offset += line_bytes;
                continue;
            }
            
            if (line[0] == '>' || line[0] == '@') {
                // New sequence
                if (seq_len > 0) {
                    // Save the previous sequence
                    if (format_ == FileFormat::FASTA) {
                        output << name << "\t" << seq_len << "\t" << seq_offset 
                              << "\t" << line_bases << "\t" << line_len << "\n";
                    } else {
                        output << name << "\t" << seq_len << "\t" << seq_offset 
                              << "\t" << line_bases << "\t" << line_len 
                              << "\t" << qual_offset << "\n";
                    }
                    
                    // Reset for new sequence
                    seq_len = 0;
                    line_len = 0;
                    line_bases = 0;
                }
                
                // Extract name from header line
                name = line.substr(1, line.find_first_of(" \t") - 1);
                
                // Set format based on first character
                if (line[0] == '>') {
                    format_ = FileFormat::FASTA;
                    state = IN_SEQ;
                } else {
                    format_ = FileFormat::FASTQ;
                    state = IN_SEQ;
                }
                
                // Remember offset of sequence
                seq_offset = offset + line_bytes;
            } else if (line[0] == '+' && state == IN_SEQ && format_ == FileFormat::FASTQ) {
                // Quality header
                state = IN_QUAL;
                qual_offset = offset + line_bytes;
            } else if (state == IN_SEQ) {
                // Sequence line
                if (line_len == 0) {
                    // First line of sequence
                    line_len = line_bytes;
                    line_bases = line.length();
                } else if (line_bytes != line_len) {
                    // Inconsistent line length
                    return false;
                }
                
                seq_len += line.length();
            } else if (state == IN_QUAL) {
                // Quality line
                if (line_bytes != line_len) {
                    // Inconsistent line length
                    return false;
                }
            }
            
            offset += line_bytes;
        }
        
        // Save the last sequence
        if (seq_len > 0) {
            if (format_ == FileFormat::FASTA) {
                output << name << "\t" << seq_len << "\t" << seq_offset 
                      << "\t" << line_bases << "\t" << line_len << "\n";
            } else {
                output << name << "\t" << seq_len << "\t" << seq_offset 
                      << "\t" << line_bases << "\t" << line_len 
                      << "\t" << qual_offset << "\n";
            }
        }
        
        return true;
    }
    
public:
    /**
     * @brief Construct a new FastaReader
     * 
     * @param fasta_path Path to the FASTA/FASTQ file
     * @param build_index Whether to build an index if it doesn't exist
     */
    FastaReader(const std::string& fasta_path, bool build_index = false) 
        : filename_(fasta_path) {
        
        // Check if file exists first
        FILE* test_file = fopen(fasta_path.c_str(), "r");
        if (!test_file) {
            throw std::runtime_error("Cannot open file: " + fasta_path + " - " + strerror(errno));
        }
        fclose(test_file);
        
        try {
            std::cout << "Opening FASTA file: " << fasta_path << std::endl;
            
            int flags = 0;
            if (build_index) {
                flags |= FAI_CREATE;
            }
            
            // Try to use htslib faidx directly
            faidx_t* fai = fai_load3(fasta_path.c_str(), NULL, NULL, flags);
            
            if (!fai) {
                std::string error_msg = "Failed to load or build index for " + fasta_path;
                if (errno != 0) {
                    error_msg += " - " + std::string(strerror(errno));
                }
                throw std::runtime_error(error_msg);
            }
            
            // Successfully loaded with htslib, extract the index data
            format_ = FileFormat::FASTA; // Default to FASTA, we'll check for FASTQ later
            
            // Get sequence names and lengths
            int n_seqs = faidx_nseq(fai);
            std::cout << "Found " << n_seqs << " sequences in index" << std::endl;
            
            for (int i = 0; i < n_seqs; i++) {
                const char* name = faidx_iseq(fai, i);
                sequence_names_.push_back(name);
                
                // Create an entry
                IndexEntry entry;
                entry.name = name;
                entry.length = faidx_seq_len64(fai, name);
                
                // Open the index file to get more details
                std::string fai_path = fasta_path + ".fai";
                std::ifstream index_file(fai_path);
                if (index_file.is_open()) {
                    std::string line;
                    while (std::getline(index_file, line)) {
                        std::istringstream iss(line);
                        std::string seq_name;
                        int64_t len, offset;
                        int32_t line_bases, line_width;
                        
                        if (!(iss >> seq_name >> len >> offset >> line_bases >> line_width)) {
                            continue;  // Skip malformed lines
                        }
                        
                        if (seq_name == name) {
                            entry.offset = offset;
                            entry.line_bases = line_bases;
                            entry.line_width = line_width;
                            break;
                        }
                    }
                }
                
                entries_[name] = entry;
            }
            
            std::cout << "Index loaded successfully" << std::endl;
            
            // Remember to clean up
            fai_destroy(fai);
            
        } catch (const std::exception& e) {
            throw std::runtime_error("Error initializing FastaReader: " + std::string(e.what()));
        }
    }
    
    /**
     * @brief Load an existing index from file
     * 
     * @param fai_path Path to the .fai index file
     * @return true if successful
     */
    bool load_index(const std::string& fai_path) {
        std::ifstream index(fai_path);
        if (!index.is_open()) {
            return false;
        }
        
        // Clear existing data
        entries_.clear();
        sequence_names_.clear();
        
        // Parse each line
        std::string line;
        while (std::getline(index, line)) {
            std::istringstream iss(line);
            std::string name;
            IndexEntry entry;
            
            // Read fields
            if (!(iss >> name >> entry.length >> entry.offset >> entry.line_bases >> entry.line_width)) {
                continue; // Skip malformed lines
            }
            
            // Check if it's a FASTQ index
            if (iss >> entry.qual_offset) {
                format_ = FileFormat::FASTQ;
            } else {
                format_ = FileFormat::FASTA;
                entry.qual_offset = 0;
            }
            
            // Store the entry
            entry.name = name;
            entries_[name] = entry;
            sequence_names_.push_back(name);
        }
        
        return !entries_.empty();
    }
    
    /**
     * @brief Build a new index for the FASTA/FASTQ file
     * 
     * @param out_path Optional output path for the index file
     * @return true if successful
     */
    bool build_index(const std::string& out_path = "") {
        std::string fai_path = out_path.empty() ? filename_ + ".fai" : out_path;
        return build_index_internal(fai_path) && load_index(fai_path);
    }
    
    /**
     * @brief Retrieve sequence based on coordinates
     * 
     * @param contig Contig/chromosome name
     * @param start Start position (0-based)
     * @param end End position (0-based, exclusive)
     * @return std::string The requested sequence
     */
    std::string fetch_sequence(const std::string& contig, int64_t start, int64_t end) const {
        try {
            // Use htslib directly for more robust sequence fetching
            faidx_t* fai = fai_load(filename_.c_str());
            if (!fai) {
                throw std::runtime_error("Failed to load index for " + filename_);
            }
            
            // Convert to 1-based, inclusive coordinates for faidx_fetch_seq64
            hts_pos_t len;
            char* seq = faidx_fetch_seq64(fai, contig.c_str(), start, end-1, &len);
            
            if (!seq) {
                if (!faidx_has_seq(fai, contig.c_str())) {
                    fai_destroy(fai);
                    throw std::runtime_error("Sequence not found: " + contig);
                } else {
                    fai_destroy(fai);
                    throw std::runtime_error("Failed to fetch sequence region: " + 
                                            contig + ":" + std::to_string(start) + 
                                            "-" + std::to_string(end));
                }
            }
            
            std::string result(seq, len);
            free(seq);
            fai_destroy(fai);
            
            return result;
        } catch (const std::exception& e) {
            throw std::runtime_error("Error fetching sequence: " + std::string(e.what()));
        }
    }
    
    /**
     * @brief Retrieve sequence based on region string
     * 
     * @param region Region string (e.g., "chr1:1000-2000")
     * @return std::string The requested sequence
     */
    std::string fetch_sequence(const std::string& region) const {
        std::string contig;
        int64_t start, end;
        
        if (!parse_region(region, contig, start, end)) {
            throw std::runtime_error("Invalid region format: " + region);
        }
        
        return fetch_sequence(contig, start, end);
    }
    
    /**
     * @brief Retrieve quality scores based on coordinates (FASTQ only)
     * 
     * @param contig Contig/chromosome name
     * @param start Start position (0-based)
     * @param end End position (0-based, exclusive)
     * @return std::string The requested quality scores
     */
    std::string fetch_quality(const std::string& contig, int64_t start, int64_t end) const {
        if (format_ != FileFormat::FASTQ) {
            throw std::runtime_error("Quality scores only available for FASTQ files");
        }
        
        // Find the entry
        auto it = entries_.find(contig);
        if (it == entries_.end()) {
            throw std::runtime_error("Sequence not found: " + contig);
        }
        
        const IndexEntry& entry = it->second;
        
        // Boundary checks
        if (start < 0) start = 0;
        if (end > entry.length) end = entry.length;
        if (start >= end) return "";
        
        // Retrieve the quality scores
        return retrieve_data(entry, entry.qual_offset, start, end);
    }
    
    /**
     * @brief Retrieve quality scores based on region string (FASTQ only)
     * 
     * @param region Region string (e.g., "chr1:1000-2000")
     * @return std::string The requested quality scores
     */
    std::string fetch_quality(const std::string& region) const {
        if (format_ != FileFormat::FASTQ) {
            throw std::runtime_error("Quality scores only available for FASTQ files");
        }
        
        std::string contig;
        int64_t start, end;
        
        if (!parse_region(region, contig, start, end)) {
            throw std::runtime_error("Invalid region format: " + region);
        }
        
        return fetch_quality(contig, start, end);
    }
    
    /**
     * @brief Get all available sequence names
     * 
     * @return std::vector<std::string> List of sequence names
     */
    std::vector<std::string> get_sequence_names() const {
        return sequence_names_;
    }
    
    /**
     * @brief Get sequence length for a contig
     * 
     * @param contig Contig/chromosome name
     * @return int64_t Length of the sequence
     */
    int64_t get_sequence_length(const std::string& contig) const {
        auto it = entries_.find(contig);
        if (it == entries_.end()) {
            return -1; // Sequence not found
        }
        return it->second.length;
    }
    
    /**
     * @brief Get the file format
     * 
     * @return FileFormat Format of the file
     */
    FileFormat get_format() const {
        return format_;
    }
    
    /**
     * @brief Destructor handles cleanup
     */
    ~FastaReader() {
        // Close all thread-local file handles
        std::lock_guard<std::mutex> lock(resource_mutex_);
        thread_local_files_.clear();
    }
};

} // namespace ts_faidx

#endif // THREAD_SAFE_FAIDX_HPP
