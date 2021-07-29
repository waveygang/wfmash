/**
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

#include "map/include/map_parameters.hpp"
#include "map/include/base_types.hpp"
#include "map/include/winSketch.hpp"
#include "map/include/computeMap.hpp"
#include "map/include/parseCmdArgs.hpp"

#include "align/include/align_parameters.hpp"
#include "align/include/computeAlignments.hpp"
#include "align/include/parseCmdArgs.hpp"

#include "yeet/include/parse_args.hpp"

//External includes
#include "common/args.hxx"
#include "common/ALeS.hpp"

int main(int argc, char** argv) {
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
        skch::printCmdOptions(map_parameters);

        auto t0 = skch::Time::now();

        if (map_parameters.use_spaced_seeds) {
          std::cerr << "[wfmash::map] Generating spaced seeds" << std::endl;
          uint32_t seed_weight = map_parameters.spaced_seed_params.weight;
          uint32_t seed_count = map_parameters.spaced_seed_params.seed_count;
          float similarity = map_parameters.spaced_seed_params.similarity;
          uint32_t region_length = map_parameters.spaced_seed_params.region_length;

          ales::spaced_seeds sps = ales::generate_spaced_seeds(seed_weight, seed_count, similarity, region_length);
          std::chrono::duration<double> time_spaced_seeds = skch::Time::now() - t0;
          std::cerr << "[wfmash::map] Time spent generating spaced seeds " << time_spaced_seeds.count()  << " seconds" << std::endl;
          map_parameters.spaced_seed_sensitivity = sps.sensitivity;
          map_parameters.spaced_seeds =  sps.seeds;
          ales::printSpacedSeeds(map_parameters.spaced_seeds);
          std::cerr << "[wfmash::map] Spaced seed sensitivity " << sps.sensitivity << std::endl;
        }

        //Build the sketch for reference
        skch::Sketch referSketch(map_parameters);

        std::chrono::duration<double> timeRefSketch = skch::Time::now() - t0;
        std::cerr << "[wfmash::map] time spent computing the reference index: " << timeRefSketch.count() << " sec" << std::endl;

        //Map the sequences in query file
        t0 = skch::Time::now();

        skch::Map mapper = skch::Map(map_parameters, referSketch);

        std::chrono::duration<double> timeMapQuery = skch::Time::now() - t0;
        std::cerr << "[wfmash::map] time spent mapping the query: " << timeMapQuery.count() << " sec" << std::endl;
        std::cerr << "[wfmash::map] mapping results saved in: " << map_parameters.outFileName << std::endl;

        if (yeet_parameters.approx_mapping) {
            return 0;
        }

        if (align_parameters.sam_format) {
            std::ofstream outstrm(align_parameters.pafOutputFile);
            for (auto& x : referSketch.metadata){
                outstrm << "@SQ\tSN:" << x.name << "\tLN:" << x.len << "\n";
            }
            outstrm << "@PG\tID:wfmash\tPN:wfmash\tVN:0.1\tCL:wfmash\n";
            outstrm.close();
        }
     } else {
        robin_hood::unordered_flat_map< std::string, skch::seqno_t > seqName_to_seqCounter;
        robin_hood::unordered_flat_map< std::string, uint64_t > seqName_to_seqLen;

        //sequence counter while parsing file
        skch::seqno_t seqCounter = 0;

        for(const auto &fileName : map_parameters.refSequences)
        {
            seqiter::for_each_seq_in_file(
                    fileName,
                    [&](const std::string& seq_name,
                            const std::string& seq) {
                        seqName_to_seqCounter[seq_name] = seqCounter++;
                        seqName_to_seqLen[seq_name] = seq.length();
                    });
        }

        std::ifstream mappingListStream(map_parameters.outFileName);
        std::string mappingRecordLine;
        align::MappingBoundaryRow currentRecord;
        std::vector<align::MappingBoundaryRow> allReadMappings;  //Aggregate mapping results for the complete run

        while (!mappingListStream.eof()){
            std::getline(mappingListStream, mappingRecordLine);
            if( !mappingRecordLine.empty() ) {
                parseMashmapRow(mappingRecordLine, currentRecord);

                allReadMappings.push_back(currentRecord);
            }
        }

        std::sort(allReadMappings.begin(), allReadMappings.end(), [&seqName_to_seqCounter](const align::MappingBoundaryRow &a, const align::MappingBoundaryRow &b)
        {
            return (seqName_to_seqCounter[a.qId] < seqName_to_seqCounter[b.qId]);
        });

        std::ofstream outstrm(align_parameters.mashmapPafFile);
        for(auto &e : allReadMappings)
        {
            outstrm << e.qId
            << "\t" << seqName_to_seqLen[e.qId]
            << "\t" << e.qStartPos
            << "\t" << e.qEndPos
            << "\t" << (e.strand == skch::strnd::FWD ? "+" : "-")
            << "\t" << e.refId
            << "\t" << seqName_to_seqLen[e.refId]
            << "\t" << e.rStartPos
            << "\t" << e.rEndPos
            << "\t" << "0"
            << "\t" << "0"
            << "\t" << "0"
            << "\t" << "id:f:" << e.mashmap_identity * 100.0
            << "\n";
        }
    }

    align::printCmdOptions(align_parameters);
    align::Aligner alignObj(align_parameters);

    auto t0 = skch::Time::now();
    std::chrono::duration<double> timeRefRead = skch::Time::now() - t0;
    std::cerr << "[wfmash::align] time spent read the reference sequences: " << timeRefRead.count() << " sec" << std::endl;

    //Compute the alignments
    alignObj.compute();

    std::chrono::duration<double> timeAlign = skch::Time::now() - t0;
    std::cerr << "[wfmash::align] time spent computing the alignment: " << timeAlign.count() << " sec" << std::endl;

    std::cerr << "[wfmash::align] alignment results saved in: " << align_parameters.pafOutputFile << std::endl;

}
