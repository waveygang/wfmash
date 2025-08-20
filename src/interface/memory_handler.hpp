#ifndef WFMASH_MEMORY_HANDLER_HPP
#define WFMASH_MEMORY_HANDLER_HPP

#include <new>
#include <iostream>
#include <atomic>
#include <chrono>
#include <thread>
#include <sstream>

namespace wfmash {
namespace memory {

// Track allocation failures for informative messages
inline std::atomic<int> allocation_failures(0);
inline std::atomic<bool> handler_installed(false);
inline std::chrono::steady_clock::time_point last_failure_time;

// Store actual parameters for better error messages
inline int actual_threads = 0;
inline std::string actual_batch_size = "";

inline void memory_exhausted_handler() {
    // Print error message immediately (no rate limiting since we're exiting)
    std::cerr << "\n========================================\n";
    std::cerr << "[wfmash] ERROR: Memory allocation failed!\n";
    std::cerr << "========================================\n\n";
    
    std::cerr << "The system has run out of available memory.\n";
    std::cerr << "This typically happens when:\n";
    std::cerr << "  - The batch size (-b) is too large for available RAM\n";
    std::cerr << "  - Too many threads are trying to allocate memory simultaneously\n";
    std::cerr << "  - The input sequences are very large\n\n";
    
    std::cerr << "SUGGESTIONS TO FIX THIS:\n";
    
    // Provide specific suggestions based on current settings
    if (!actual_batch_size.empty()) {
        std::cerr << "  1. Reduce batch size: Current is -b " << actual_batch_size << "\n";
        std::cerr << "     Try: -b 500m, -b 100m, or -b 50m\n";
    } else {
        std::cerr << "  1. Reduce batch size: Try -b 500m or -b 100m instead of -b 1g\n";
    }
    
    if (actual_threads > 0) {
        std::cerr << "  2. Use fewer threads: Current is -t " << actual_threads << "\n";
        if (actual_threads > 48) {
            std::cerr << "     Try: -t " << (actual_threads / 2) << " or -t 24\n";
        } else if (actual_threads > 8) {
            std::cerr << "     Try: -t " << (actual_threads / 2) << " or -t 8\n";
        } else {
            std::cerr << "     Try: -t " << std::max(1, actual_threads / 2) << "\n";
        }
    } else {
        std::cerr << "  2. Use fewer threads: Try -t 24 or -t 48 instead of -t 112\n";
    }
    std::cerr << "  3. Run on a node with more memory\n";
    std::cerr << "  4. Split your input into smaller chunks\n\n";
    
    std::cerr << "The program will now exit to prevent system instability.\n";
    std::cerr << "Please adjust parameters and try again.\n";
    std::cerr << "========================================\n\n";
        
    
    // Flush output
    std::cerr.flush();
    
    // Exit immediately with error code 1
    std::exit(1);
}

inline void install_memory_handler(int threads = 0, const std::string& batch_size = "") {
    if (!handler_installed.exchange(true)) {
        actual_threads = threads;
        actual_batch_size = batch_size;
        std::set_new_handler(memory_exhausted_handler);
        std::cerr << "[wfmash] Memory handler installed - will abort cleanly on allocation failure\n";
    }
}

inline void uninstall_memory_handler() {
    if (handler_installed.exchange(false)) {
        std::set_new_handler(nullptr);
    }
}

} // namespace memory
} // namespace wfmash

#endif // WFMASH_MEMORY_HANDLER_HPP