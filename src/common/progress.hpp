#pragma once

#include <iostream>
#include <string>
#include <atomic>
#include <thread>
#include <chrono>
#include <iomanip>

namespace progress_meter {

class ProgressMeter {
private:
    const uint64_t update_interval = 500; // ms between updates
    const uint64_t min_progress_for_update = 1000; // Minimum progress before showing an update
    std::atomic<bool> running;
    std::chrono::time_point<std::chrono::steady_clock> last_update;

public:
    std::string banner;
    std::atomic<uint64_t> total;
    std::atomic<uint64_t> completed;
    std::chrono::time_point<std::chrono::steady_clock> start_time;
    std::thread logger;
    ProgressMeter(uint64_t _total, const std::string& _banner)
        : total(_total), banner(_banner), running(true) {
        start_time = std::chrono::steady_clock::now();
        last_update = start_time;
        completed = 0;
        
        logger = std::thread([this]() {
            uint64_t last_completed = 0;
            
            while (running.load(std::memory_order_relaxed)) {
                auto now = std::chrono::steady_clock::now();
                auto time_since_update = std::chrono::duration_cast<std::chrono::milliseconds>(now - last_update).count();
                uint64_t current_completed = completed.load(std::memory_order_relaxed);
                
                if (time_since_update >= update_interval && 
                    (current_completed - last_completed >= min_progress_for_update ||
                     current_completed >= total)) {
                    do_print();
                    last_completed = current_completed;
                    last_update = now;
                }
                
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
            }
        });
    };
    void do_print(void) {
        auto curr = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = curr-start_time;
        double rate = completed / elapsed_seconds.count();
        double seconds_to_completion = 0;
        if (completed > 0) {
            if (completed >= total) {
                seconds_to_completion = 0;
            } else {
                seconds_to_completion = (total - completed) / rate;
                // Cap unreasonable estimates at 100 days
                if (seconds_to_completion > 8640000) { // 100 days in seconds
                    seconds_to_completion = 8640000;
                }
                // Prevent negative estimates
                if (seconds_to_completion < 0) {
                    seconds_to_completion = 0;
                }
            }
        }
        std::cerr << "\r" << banner << " "
                  << std::defaultfloat
                  << std::setfill(' ')
                  << std::setw(5)
                  << std::fixed
                  << std::setprecision(2)
                  << 100.0 * ((double)completed / (double)total) << "% "
                  << "in: " << print_time(elapsed_seconds.count()) << " "
                  << "todo: " << print_time(seconds_to_completion) << " @"
                  << std::setw(4) << std::scientific << rate << "/s";
    }
    void finish() {
        running.store(false, std::memory_order_relaxed);
        if (logger.joinable()) {
            logger.join();
        }
        completed.store(total);
        do_print();
        std::cerr << std::endl;
    }
    std::string print_time(const double& _seconds) {
        int days = 0, hours = 0, minutes = 0, seconds = 0;
        distribute_seconds(days, hours, minutes, seconds, _seconds);
        std::stringstream buffer;
        buffer << std::setfill('0') << std::setw(2) << days << ":"
               << std::setfill('0') << std::setw(2) << hours << ":"
               << std::setfill('0') << std::setw(2) << minutes << ":"
               << std::setfill('0') << std::setw(2) << seconds;
        return buffer.str();
    }
    void distribute_seconds(int& days, int& hours, int& minutes, int& seconds, const double& input_seconds) {
        const int cseconds_in_day = 86400;
        const int cseconds_in_hour = 3600;
        const int cseconds_in_minute = 60;
        const int cseconds = 1;
        days = std::floor(input_seconds / cseconds_in_day);
        hours = std::floor(((int)input_seconds % cseconds_in_day) / cseconds_in_hour);
        minutes = std::floor((((int)input_seconds % cseconds_in_day) % cseconds_in_hour) / cseconds_in_minute);
        seconds = ((((int)input_seconds % cseconds_in_day) % cseconds_in_hour) % cseconds_in_minute) / cseconds; // + (input_seconds - std::floor(input_seconds));
        //std::cerr << input_seconds << " seconds is " << days << " days, " << hours << " hours, " << minutes << " minutes, and " << seconds << " seconds." << std::endl;
    }
    void increment(const uint64_t& incr) {
        completed.fetch_add(incr, std::memory_order_relaxed);
    }
    ~ProgressMeter() {
        if (running.load(std::memory_order_relaxed)) {
            finish();
        }
    }
};

}
