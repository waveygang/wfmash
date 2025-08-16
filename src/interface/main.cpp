/**
 * This is the main entry point for wfmash
 *
 * @file    align.cpp
 * @ingroup src
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#include <sstream>
#include <fstream>
#include <iostream>
#include <ctime>
#include <chrono>
#include <functional>
#include <cstdio>
#include <iomanip>
#include <memory>
#ifdef __GLIBC__
#include <malloc.h>
#endif

#include "map/include/computeMap.hpp"
#include "map/include/parseCmdArgs.hpp"
#include "map/include/sequenceIds.hpp"
#include "map/include/map_stats.hpp"

#include "map/include/winSketch.hpp"

#include "interface/parse_args.hpp"
#include "align/include/parseCmdArgs.hpp"
#include "align/include/computeAlignments.hpp"


// External includes
#include "common/args.hxx"
#include "common/ALeS.hpp"

// Memory allocation failure handler
#include <new>
#include <atomic>
#include <thread>

void wfmash_memory_handler() {
    static std::atomic<int> attempts(0);
    int attempt_count = ++attempts;
    
    if (attempt_count > 3) {
        std::cerr << "\n[wfmash::fatal] Memory allocation failed after " 
                  << attempt_count << " attempts\n";
        std::cerr << "The system has run out of memory or reached allocation limits.\n";
        std::cerr << "\nSuggestions to resolve:\n";
        std::cerr << "  1. Use -b flag with smaller value (e.g., -b 50m or -b 100m)\n";
        std::cerr << "  2. Reduce thread count with -t (fewer concurrent tasks)\n";
        std::cerr << "  3. Process smaller sequence subsets\n";
        std::cerr << "  4. Run on a system with more available memory\n";
        std::cerr << "\nAborting to prevent incomplete results.\n";
        std::abort();
    }
    
    // Log warning for first few attempts
    std::cerr << "\n[wfmash::warning] Memory allocation failed (attempt " 
              << attempt_count << "/3)\n";
    std::cerr << "Retrying allocation after brief pause...\n";
    
    // Brief pause to allow system to free memory
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
}

int main(int argc, char** argv) {
    // Install memory allocation failure handler
    std::set_new_handler(wfmash_memory_handler);
    /*
     * Make sure env variable MALLOC_ARENA_MAX is unset
     * for efficient multi-thread execution
     */
    unsetenv((char *)"MALLOC_ARENA_MAX");

    // get our parameters from the command line
    skch::Parameters map_parameters;
    align::Parameters align_parameters;
    yeet::Parameters yeet_parameters;
    yeet::parse_args(argc, argv, map_parameters, align_parameters, yeet_parameters);

    //parameters.refSequences.push_back(ref);

    //skch::parseandSave(argc, argv, cmd, parameters);
    if (!yeet_parameters.remapping) {
        // Handle auto percentage identity estimation
        if (map_parameters.auto_pct_identity) {
            std::cerr << "[wfmash] ANI-based identity estimation enabled (ani" << map_parameters.ani_percentile;
            if (map_parameters.ani_adjustment != 0) {
                std::cerr << std::showpos << map_parameters.ani_adjustment << std::noshowpos;
            }
            std::cerr << ")..." << std::endl;

            // Instantiate the SequenceIdManager here to get the correct grouping information.
            // This is the single source of truth for grouping.
            // Create target prefix vector, handling empty strings properly
            std::vector<std::string> target_prefix_vec;
            if (!map_parameters.target_prefix.empty()) {
                target_prefix_vec.push_back(map_parameters.target_prefix);
            }
            
            auto idManager = std::make_unique<skch::SequenceIdManager>(
                map_parameters.querySequences,
                map_parameters.refSequences,
                map_parameters.query_prefix,
                target_prefix_vec,
                std::string(1, map_parameters.prefix_delim),
                map_parameters.query_list,
                map_parameters.target_list
            );

            try {
                // Call the estimation function, which will respect the groups defined in the idManager.
                double estimated_identity = skch::Stat::estimate_identity_for_groups(map_parameters, *idManager);
                
                // Update the parameters that the rest of the program will use.
                map_parameters.percentageIdentity = estimated_identity;
                
                // Recalculate sketch size based on the new identity threshold
                // Only if sketch size was auto-calculated (not manually specified with -s)
                if (!map_parameters.sketch_size_manually_set) {
                    const double md = 1 - map_parameters.percentageIdentity;
                    double dens = 0.02 * (1 + (md / 0.1));
                    int old_sketch_size = map_parameters.sketchSize;
                    map_parameters.sketchSize = dens * (map_parameters.windowLength - map_parameters.kmerSize);
                    
                    // Ensure sketch size doesn't exceed window size
                    if (map_parameters.sketchSize > map_parameters.windowLength) {
                        map_parameters.sketchSize = map_parameters.windowLength;
                    }
                    
                    std::cerr << "[wfmash] Updated sketch size from " << old_sketch_size 
                              << " to " << map_parameters.sketchSize 
                              << " based on estimated identity" << std::endl;
                }
                
                std::cerr << "[wfmash] Using estimated identity cutoff: " 
                          << std::fixed << std::setprecision(2) << estimated_identity * 100 << "%" << std::endl;
            } catch (const std::exception& e) {
                std::cerr << "[wfmash] Error during identity estimation: " << e.what() << std::endl;
                std::cerr << "[wfmash] Falling back to default identity threshold: " 
                          << std::fixed << std::setprecision(2) << skch::fixed::percentage_identity * 100 << "%" << std::endl;
                map_parameters.percentageIdentity = skch::fixed::percentage_identity;
                map_parameters.auto_pct_identity = false;
            }
        }

        skch::printCmdOptions(map_parameters);
        
        // Log final parameters after ANI estimation
        int minimum_hits = std::max(map_parameters.minimum_hits, 
                                   skch::Stat::estimateMinimumHitsRelaxed(map_parameters.sketchSize, 
                                                                          map_parameters.kmerSize, 
                                                                          map_parameters.percentageIdentity, 
                                                                          skch::fixed::confidence_interval));
        
        std::cerr << "[wfmash] Final parameters: identity=" << std::fixed << std::setprecision(1) 
                  << map_parameters.percentageIdentity * 100 << "%, "
                  << "sketchSize=" << map_parameters.sketchSize << ", "
                  << "minimumHits=" << minimum_hits << ", "
                  << "windowLength=" << map_parameters.windowLength << ", "
                  << "kmerSize=" << map_parameters.kmerSize << std::endl;

        auto t0 = skch::Time::now();

        if (map_parameters.use_spaced_seeds) {
          std::cerr << "[wfmash::mashmap] Generating spaced seeds..." << std::endl;
          uint32_t seed_weight = map_parameters.spaced_seed_params.weight;
          uint32_t seed_count = map_parameters.spaced_seed_params.seed_count;
          float similarity = map_parameters.spaced_seed_params.similarity;
          uint32_t region_length = map_parameters.spaced_seed_params.region_length;

          ales::spaced_seeds sps = ales::generate_spaced_seeds(seed_weight, seed_count, similarity, region_length);
          std::chrono::duration<double> time_spaced_seeds = skch::Time::now() - t0;
          map_parameters.spaced_seed_sensitivity = sps.sensitivity;
          map_parameters.spaced_seeds =  sps.seeds;
          ales::printSpacedSeeds(map_parameters.spaced_seeds);
          std::cerr << "[wfmash::mashmap] Generated spaced seeds in " << time_spaced_seeds.count() << "s (sensitivity: " << sps.sensitivity << ")" << std::endl;
        }

        //Map the sequences in query file
        t0 = skch::Time::now();

        skch::Map mapper = skch::Map(map_parameters);

        std::chrono::duration<double> timeMapQuery = skch::Time::now() - t0;
        std::cerr << "[wfmash::mashmap] Mapped query in " << timeMapQuery.count() << "s, results saved to: " << map_parameters.outFileName << std::endl;

        if (yeet_parameters.approx_mapping) {
            return 0;
        }
        
        // Trim memory after mapping phase to release unused memory back to OS
        #ifdef __GLIBC__
        malloc_trim(0);
        #endif
     }

    align::printCmdOptions(align_parameters);

    auto t0 = skch::Time::now();
    align::Aligner alignObj(align_parameters);
    std::chrono::duration<double> timeRefRead = skch::Time::now() - t0;
    std::cerr << "[wfmash::align] time spent loading the reference index: " << timeRefRead.count() << " sec" << std::endl;

    //Compute the alignments
    alignObj.compute();

    std::chrono::duration<double> timeAlign = skch::Time::now() - t0;
    std::cerr << "[wfmash::align] time spent computing the alignment: " << timeAlign.count() << " sec" << std::endl;

    std::cerr << "[wfmash::align] alignment results saved in: " << align_parameters.pafOutputFile << std::endl;

}
