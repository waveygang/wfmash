#ifndef WFMASH_MEMORY_MONITOR_HPP
#define WFMASH_MEMORY_MONITOR_HPP

#include <atomic>
#include <chrono>

namespace wfmash {
namespace memory {

// Global atomic counters for memory allocation tracking
inline std::atomic<int> tasks_stalled(0);      // Tasks waiting for memory
inline std::atomic<int> tasks_executing(0);    // Tasks currently executing
inline std::atomic<int> total_stall_events(0); // Total stall events seen

// Track when we last logged to avoid spam
inline std::atomic<std::chrono::steady_clock::time_point> last_log_time{std::chrono::steady_clock::now()};

// Helper to get memory stats string
inline std::string get_memory_status() {
    int threads_with_stalls = tasks_stalled.load();
    int total_events = total_stall_events.load();
    
    // Only show status if there have been stalls
    if (threads_with_stalls > 0 || total_events > 0) {
        return " [" + std::to_string(threads_with_stalls) + " threads experienced stalls]";
    }
    return "";
}

} // namespace memory
} // namespace wfmash

#endif // WFMASH_MEMORY_MONITOR_HPP