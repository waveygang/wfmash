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

inline void memory_exhausted_handler() {
    allocation_failures.fetch_add(1);
    auto now = std::chrono::steady_clock::now();
    
    // Only print message once per second to avoid spam
    static std::chrono::steady_clock::time_point last_print = std::chrono::steady_clock::time_point::min();
    if (std::chrono::duration_cast<std::chrono::seconds>(now - last_print).count() >= 1) {
        std::cerr << "\n========================================\n";
        std::cerr << "[wfmash] ERROR: Memory allocation failed!\n";
        std::cerr << "========================================\n\n";
        
        std::cerr << "The system has run out of available memory.\n";
        std::cerr << "This typically happens when:\n";
        std::cerr << "  - The batch size (-b) is too large for available RAM\n";
        std::cerr << "  - Too many threads are trying to allocate memory simultaneously\n";
        std::cerr << "  - The input sequences are very large\n\n";
        
        std::cerr << "SUGGESTIONS TO FIX THIS:\n";
        std::cerr << "  1. Reduce batch size: Try -b 500m or -b 100m instead of -b 1g\n";
        std::cerr << "  2. Use fewer threads: Try -t 24 or -t 48 instead of -t 112\n";
        std::cerr << "  3. Run on a node with more memory\n";
        std::cerr << "  4. Split your input into smaller chunks\n\n";
        
        std::cerr << "The program will now exit to prevent system instability.\n";
        std::cerr << "Please adjust parameters and try again.\n";
        std::cerr << "========================================\n\n";
        
        last_print = now;
    }
    
    // Give a moment for the message to be written
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    
    // Exit cleanly with error code
    std::abort();
}

inline void install_memory_handler() {
    if (!handler_installed.exchange(true)) {
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