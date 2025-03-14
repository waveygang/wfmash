#pragma once

#include <iostream>
#include <string>
#include <chrono>
#include <atomic>
#include <mutex>
#include <memory>
#include <thread>
#include <unistd.h>
#include "indicators.hpp"

namespace progress_meter {

class ProgressMeter {
private:
    std::string banner;
    std::atomic<uint64_t> total;
    std::atomic<uint64_t> completed;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
    bool use_progress_bar;
    std::atomic<bool> is_finished;
    std::atomic<bool> running;
    std::mutex mutex;
    std::unique_ptr<indicators::BlockProgressBar> progress_bar;
    std::thread update_thread;
    
    // Update interval for the progress bar (milliseconds)
    const uint64_t update_interval = 100; 

    void update_progress_thread() {
        uint64_t last_progress = 0;
        
        while (running.load()) {
            // Update progress bar at regular intervals
            if (use_progress_bar && progress_bar) {
                auto curr_progress = std::min(completed.load(), total.load());
                float progress_percent = static_cast<float>(curr_progress) / total.load() * 100.0f;
                
                // Only update if there's a meaningful change in percentage
                // or we're just starting or finishing
                if (curr_progress != last_progress && 
                    (last_progress == 0 || curr_progress == total.load() || 
                     std::abs(static_cast<float>(curr_progress - last_progress)) / total.load() >= 0.01)) {
                    
                    std::lock_guard<std::mutex> lock(mutex);
                    progress_bar->set_progress(curr_progress);
                    last_progress = curr_progress;
                    
                    // If we've reached 100%, mark as completed and stop the thread
                    if (curr_progress >= total.load()) {
                        if (!is_finished.load()) {
                            progress_bar->mark_as_completed();
                            is_finished.store(true);
                        }
                        break;  // Exit the update thread
                    }
                }
            }
            
            // Sleep for update interval
            std::this_thread::sleep_for(std::chrono::milliseconds(update_interval));
        }
    }

public:
    ProgressMeter(uint64_t _total, const std::string& _banner)
        : banner(_banner), total(_total), completed(0), is_finished(false), running(true) {
        
        start_time = std::chrono::high_resolution_clock::now();
        
        // Check if stderr is a TTY
        use_progress_bar = isatty(fileno(stderr));
        
        // Only print the banner if we're not using a progress bar
        if (!use_progress_bar) {
            std::cerr << banner << std::endl;
        }
        
        if (use_progress_bar) {
            // Hide cursor during progress display
            indicators::show_console_cursor(false);
            
            // Use the full banner as the prefix text
            // The banner should look like "[wfmash::mashmap] indexing"
            
            // Create progress bar including the full banner text
            {
                std::lock_guard<std::mutex> lock(mutex);
                progress_bar = std::make_unique<indicators::BlockProgressBar>(
                    indicators::option::BarWidth{50},
                    indicators::option::Start{"["},
                    indicators::option::End{"]"},
                    indicators::option::ShowElapsedTime{true},
                    indicators::option::ShowRemainingTime{true},
                    indicators::option::PrefixText{banner + " "},
                    // Use default color instead of blue for better visibility
                    indicators::option::MaxProgress{total},
                    indicators::option::Stream{std::cerr}
                );
            }
            
            // Start the update thread
            update_thread = std::thread(&ProgressMeter::update_progress_thread, this);
        }
    }

    ~ProgressMeter() {
        if (!is_finished.load()) {
            finish();
        }
        
        // Stop the thread
        running.store(false);
        if (update_thread.joinable()) {
            update_thread.join();
        }
        
        // Always show cursor when done
        if (use_progress_bar) {
            indicators::show_console_cursor(true);
        }
    }

    void increment(const uint64_t& incr) {
        uint64_t previous = completed.fetch_add(incr, std::memory_order_relaxed);
        uint64_t new_val = previous + incr;
        
        // For large increments, force an immediate update
        if (incr > total / 20 && use_progress_bar && progress_bar) {
            std::lock_guard<std::mutex> lock(mutex);
            progress_bar->set_progress(std::min(new_val, total.load()));
        }
    }

    void finish() {
        // Use atomic to prevent multiple finish() calls from different threads
        bool expected = false;
        if (!is_finished.compare_exchange_strong(expected, true)) {
            return; // Already finished
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        
        if (use_progress_bar && progress_bar) {
            // Final update with lock
            std::lock_guard<std::mutex> lock(mutex);
            // Directly set to 100% complete to avoid partial updates
            progress_bar->set_progress(total.load());
            progress_bar->mark_as_completed();
        }
        
        // Always print completion message
        // Only print completion message if we're not already showing a progress bar
        if (!use_progress_bar || !progress_bar) {
            std::cerr << banner << " [completed in " 
                     << elapsed.count() << "s]" << std::endl;
        }
        
        // Ensure the update thread stops
        running.store(false);
    }
};

} // namespace progress_meter
