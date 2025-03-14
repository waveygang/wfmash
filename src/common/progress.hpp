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
        while (running.load()) {
            // Update progress bar at regular intervals
            if (use_progress_bar && progress_bar) {
                std::lock_guard<std::mutex> lock(mutex);
                progress_bar->set_progress(std::min(completed.load(), total.load()));
            }
            
            // Sleep for update interval
            std::this_thread::sleep_for(std::chrono::milliseconds(update_interval));
        }
    }

public:
    ProgressMeter(uint64_t _total, const std::string& _banner)
        : banner(_banner), total(_total), completed(0), is_finished(false), running(true) {
        
        start_time = std::chrono::high_resolution_clock::now();
        
        // Always use the progress bar for debugging 
        use_progress_bar = true; // Temporarily disabled TTY check: isatty(fileno(stderr));
        
        // Print start message with timestamp
        std::time_t start_time_t = std::chrono::system_clock::to_time_t(
            std::chrono::system_clock::now());
        
        std::cerr << banner << std::endl;
                 
        if (use_progress_bar) {
            // Hide cursor during progress display
            indicators::show_console_cursor(false);
            
            // Create progress bar
            {
                std::lock_guard<std::mutex> lock(mutex);
                progress_bar = std::make_unique<indicators::BlockProgressBar>(
                    indicators::option::BarWidth{50},
                    indicators::option::Start{"["},
                    indicators::option::End{"]"},
                    indicators::option::ForegroundColor{indicators::Color::green},
                    indicators::option::ShowElapsedTime{true},
                    indicators::option::ShowRemainingTime{true},
                    indicators::option::PrefixText{banner},
                    indicators::option::FontStyles{
                        std::vector<indicators::FontStyle>{indicators::FontStyle::bold}
                    },
                    indicators::option::MaxProgress{total}
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
        completed.fetch_add(incr, std::memory_order_relaxed);
        // No need to update the bar here, the thread will handle it
    }

    void finish() {
        // Use atomic to prevent multiple finish() calls from different threads
        bool expected = false;
        if (!is_finished.compare_exchange_strong(expected, true)) {
            return; // Already finished
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        
        if (use_progress_bar) {
            // Final update with lock
            std::lock_guard<std::mutex> lock(mutex);
            if (progress_bar) {
                progress_bar->set_option(indicators::option::PrefixText{banner + " [completed]"});
                progress_bar->set_progress(total);
                progress_bar->mark_as_completed();
            }
        } else {
            // If not using progress bar, just print completion message
            std::cerr << banner << " [completed in " 
                      << elapsed.count() << "s]" << std::endl;
        }
    }
};

} // namespace progress_meter
