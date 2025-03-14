#pragma once

#include <iostream>
#include <string>
#include <chrono>
#include <atomic>
#include <mutex>
#include <memory>
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
    bool is_finished;
    std::mutex mutex;
    std::unique_ptr<indicators::BlockProgressBar> progress_bar;

public:
    ProgressMeter(uint64_t _total, const std::string& _banner)
        : banner(_banner), total(_total), completed(0), is_finished(false) {
        start_time = std::chrono::high_resolution_clock::now();
        
        // Only use progress bar if stderr is a TTY
        use_progress_bar = isatty(fileno(stderr));
        
        // Print start message with timestamp
        std::time_t start_time_t = std::chrono::system_clock::to_time_t(
            std::chrono::system_clock::now());
        
        std::cerr << banner << std::endl;
                 
        if (use_progress_bar) {
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
            
            // Hide cursor during progress display
            indicators::show_console_cursor(false);
        }
    }

    ~ProgressMeter() {
        if (!is_finished) {
            finish();
        }
        
        // Always show cursor when done
        if (use_progress_bar) {
            indicators::show_console_cursor(true);
        }
    }

    void increment(const uint64_t& incr) {
        completed.fetch_add(incr, std::memory_order_relaxed);
        if (use_progress_bar && progress_bar) {
            std::lock_guard<std::mutex> lock(mutex);
            progress_bar->set_progress(std::min(completed.load(), total.load()));
        }
    }

    void finish() {
        std::lock_guard<std::mutex> lock(mutex);
        
        if (is_finished) {
            return;
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        
        if (use_progress_bar && progress_bar) {
            progress_bar->set_option(indicators::option::PrefixText{banner + " [completed]"});
            progress_bar->set_progress(total);
            progress_bar->mark_as_completed();
            
            // Show cursor again now that we're done
            indicators::show_console_cursor(true);
        } else {
            // If not using progress bar, just print completion message
            std::cerr << banner << " [completed in " 
                      << elapsed.count() << "s]" << std::endl;
        }
        
        is_finished = true;
    }
};

} // namespace progress_meter
