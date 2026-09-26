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
#include "map/include/externalSeeder.hpp"

#include "map/include/winSketch.hpp"

#include "interface/parse_args.hpp"
#include "align/include/parseCmdArgs.hpp"
#include "align/include/computeAlignments.hpp"

// Memory handler
#include "interface/memory_handler.hpp"

// External includes
#include "common/args.hxx"
#include "common/ALeS.hpp"

// Engine handoff: the vendored low-divergence engine (the wfmash 0.14
// lineage, src/engine/low-divergence/) parses its own command line with
// the 0.14-era flags it was built for, so `--engine low-divergence ...`
// dispatches the raw remaining argv straight to it before the mainline
// parser ever sees it.
int lde_main(int argc, char** argv);

static int dispatch_engine(int argc, char** argv) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        std::string name;
        int consumed = 0;   // number of argv entries the --engine flag uses
        if (arg == "--engine" && i + 1 < argc) {
            name = argv[i + 1];
            consumed = 2;
        } else if (arg.rfind("--engine=", 0) == 0) {
            name = arg.substr(std::string("--engine=").size());
            consumed = 1;
        } else {
            continue;
        }
        if (name == "low-divergence") {
            // strip only the --engine pair; the engine parses the rest
            char** engine_argv = new char*[argc - consumed + 1];
            int j = 0;
            engine_argv[j++] = argv[0];
            for (int k = 1; k < argc; ++k) {
                if (k == i || (consumed == 2 && k == i + 1)) continue;
                engine_argv[j++] = argv[k];
            }
            int code = lde_main(j, engine_argv);
            delete[] engine_argv;
            return code;
        }
        std::cerr << "[wfmash] ERROR: unknown engine '" << name
                  << "' (available: low-divergence)" << std::endl;
        return 2;
    }
    return -1;  // no engine selection: mainline 0.24 engine
}

int main(int argc, char** argv) {
    /*
     * Make sure env variable MALLOC_ARENA_MAX is unset
     * for efficient multi-thread execution
     */
    unsetenv((char *)"MALLOC_ARENA_MAX");

    int engine_code = dispatch_engine(argc, argv);
    if (engine_code >= 0) return engine_code;

    // get our parameters from the command line
    skch::Parameters map_parameters;
    align::Parameters align_parameters;
    yeet::Parameters yeet_parameters;
    yeet::parse_args(argc, argv, map_parameters, align_parameters, yeet_parameters);
    
    // Install memory handler for clean abort on OOM
    // Convert batch size to human-readable format
    std::string batch_str;
    if (map_parameters.index_by_size == std::numeric_limits<int64_t>::max()) {
        batch_str = "unlimited";
    } else if (map_parameters.index_by_size >= 1000000000) {
        batch_str = std::to_string(map_parameters.index_by_size / 1000000000) + "g";
    } else if (map_parameters.index_by_size >= 1000000) {
        batch_str = std::to_string(map_parameters.index_by_size / 1000000) + "m";
    } else if (map_parameters.index_by_size >= 1000) {
        batch_str = std::to_string(map_parameters.index_by_size / 1000) + "k";
    } else {
        batch_str = std::to_string(map_parameters.index_by_size);
    }
    wfmash::memory::install_memory_handler(map_parameters.threads, batch_str);

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

        if (map_parameters.use_external_seeds) {
            // Process external PAF seeds instead of running MinHash mapping
            std::cerr << "[wfmash::mashmap] Using external seeds from: " << map_parameters.external_seeds_file << std::endl;

            // Create SequenceIdManager without building index (just for name->ID mapping)
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

            // Open output stream
            std::ofstream outFile;
            std::ostream* outStream = &std::cout;
            if (map_parameters.outFileName != "/dev/stdout" && map_parameters.outFileName != "-") {
                outFile.open(map_parameters.outFileName);
                if (!outFile.is_open()) {
                    std::cerr << "Error: Could not open output file: " << map_parameters.outFileName << std::endl;
                    exit(1);
                }
                outStream = &outFile;
            }

            // Process external seeds through chaining/filtering pipeline
            skch::ExternalSeeder::processExternalSeeds(
                map_parameters,
                map_parameters.external_seeds_file,
                *idManager,
                *outStream
            );

            if (outFile.is_open()) {
                outFile.close();
            }
        } else {
            // Normal MinHash-based mapping
            skch::Map mapper = skch::Map(map_parameters);
        }

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
