#pragma once

#include <iostream>
#include <string>
#include <chrono>
#include <atomic>
#include <mutex>
#include <memory>
#include <thread>
#include <unistd.h>
#include <iomanip>
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
    std::unique_ptr<indicators::BlockProgressBar> progress_bar;
    std::thread update_thread;
    
    // Update intervals (milliseconds)
    const uint64_t update_interval = 100;  // For progress bar (TTY)
    const uint64_t file_update_interval = 60000;  // 60 seconds for file output
    std::chrono::time_point<std::chrono::high_resolution_clock> last_file_update;

    void update_progress_thread() {
        uint64_t last_progress = 0;
        
        while (running.load()) {
            auto curr_time = std::chrono::high_resolution_clock::now();
            auto curr_progress = std::min(completed.load(), total.load());
            float progress_percent = static_cast<float>(curr_progress) / total.load() * 100.0f;
            
            // Update progress bar at regular intervals
            if (use_progress_bar && progress_bar) {
                // Only update if there's a meaningful change in percentage
                // or we're just starting or finishing
                if (curr_progress != last_progress && 
                    (last_progress == 0 || curr_progress == total.load() || 
                     std::abs(static_cast<float>(curr_progress - last_progress)) / total.load() >= 0.01)) {
                    
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
            } else {
                // For file output, print periodic updates or on completion
                auto elapsed_since_update = std::chrono::duration_cast<std::chrono::milliseconds>(
                    curr_time - last_file_update).count();
                
                // Update if: 1) 60 seconds have passed, 2) this is the first update, or 3) we're at 100%
                if ((elapsed_since_update >= file_update_interval || last_progress == 0 || 
                     curr_progress >= total.load()) && curr_progress != last_progress) {
                    
                    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
                        curr_time - start_time).count();
                    
                    std::cerr << banner << " [" 
                              << std::fixed << std::setprecision(1) << progress_percent << "% complete, " 
                              << curr_progress << "/" << total.load() 
                              << " units, " << elapsed << "s elapsed]" << std::endl;
                    
                    last_progress = curr_progress;
                    last_file_update = curr_time;
                    
                    // If we've reached 100%, mark as completed and exit
                    if (curr_progress >= total.load()) {
                        is_finished.store(true);
                        break;
                    }
                }
            }
            
            // If we're marked as finished, exit the loop
            if (is_finished.load()) {
                break;
            }
            
            // Sleep for update interval
            std::this_thread::sleep_for(std::chrono::milliseconds(update_interval));
        }
    }

public:
    ProgressMeter(uint64_t _total, const std::string& _banner, const bool& _use_progress_bar)
        : banner(_banner), total(_total), completed(0), is_finished(false), running(true) {
        
        start_time = std::chrono::high_resolution_clock::now();
        last_file_update = start_time;
        
        // Check if stderr is a TTY
        use_progress_bar = _use_progress_bar && isatty(fileno(stderr));
        
        // For file output, print initial banner with 0% progress
        if (!use_progress_bar) {
            std::cerr << banner << " [0.0% complete, 0/" << total.load() 
                      << " units, 0s elapsed]" << std::endl;
        }
        
        if (use_progress_bar) {
            // Hide cursor during progress display
            indicators::show_console_cursor(false);
            
            // Use the full banner as the prefix text
            // The banner should look like "[wfmash::mashmap] indexing"
            
            // Create progress bar including the full banner text
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
            
            // Start the update thread
            update_thread = std::thread(&ProgressMeter::update_progress_thread, this);
        }
    }

    ~ProgressMeter() {
        if (!is_finished.load()) {
            finish();
        } else {
            // Make sure the thread is properly stopped and joined
            running.store(false);
            if (update_thread.joinable()) {
                update_thread.join();
            }
        }
        
        // Always show cursor when done
        if (use_progress_bar) {
            indicators::show_console_cursor(true);
        }
    }

    void increment(const uint64_t& incr) {
        uint64_t previous = completed.fetch_add(incr, std::memory_order_relaxed);
        uint64_t new_val = previous + incr;
        
        // For file output with quick tasks, we might need to force an update
        // if this is a significant increment (more than 5% of total)
        if (!use_progress_bar && (incr > total.load() / 20)) {
            // Take mutex to avoid concurrent output with update thread
            static std::mutex increment_mutex;
            std::lock_guard<std::mutex> lock(increment_mutex);
            
            auto now = std::chrono::high_resolution_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start_time);
            auto elapsed_since_update = std::chrono::duration_cast<std::chrono::milliseconds>(
                now - last_file_update).count();
                
            // Only update if at least 2 seconds have passed since last update to prevent spamming
            if (elapsed_since_update > 2000) {
                float progress_percent = static_cast<float>(new_val) / total.load() * 100.0f;
                std::cerr << banner << " [" 
                          << std::fixed << std::setprecision(1) << progress_percent << "% complete, " 
                          << new_val << "/" << total.load() 
                          << " units, " << elapsed.count() << "s elapsed]" << std::endl;
                last_file_update = now;
            }
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
        
        // Set final value atomically
        completed.store(total.load(), std::memory_order_relaxed);
        
        // Force immediate update of the progress bar with the final value
        if (use_progress_bar && progress_bar) {
            progress_bar->set_progress(total.load());
            progress_bar->mark_as_completed();
        } else {
            // For file output, always print 100% completion message
            std::cerr << banner << " [100.0% complete, " << total.load() << "/" << total.load() 
                      << " units, " << elapsed.count() << "s elapsed]" << std::endl;
            std::cerr << banner << " completed." << std::endl;
        }

        // Ensure the update thread stops
        running.store(false);

        // Wait for the update thread to fully terminate before continuing
        if (update_thread.joinable()) {
            update_thread.join();
        }
        
        // Flush stderr to ensure all output is properly displayed
        std::cerr.flush();
    }
};

} // namespace progress_meter
