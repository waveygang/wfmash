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
    int stalled = tasks_stalled.load();
    int total_events = total_stall_events.load();
    
    // Only show status if there are stalls
    if (stalled > 0 || total_events > 0) {
        std::string status = " [";
        if (stalled > 0) {
            status += std::to_string(stalled) + " threads stalled";
        }
        if (total_events > 0) {
            if (stalled > 0) status += ", ";
            status += std::to_string(total_events) + " total stalls";
        }
        status += "]";
        return status;
    }
    return "";
}

} // namespace memory
} // namespace wfmash

#endif // WFMASH_MEMORY_MONITOR_HPP