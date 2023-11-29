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

#include "interface/parse_args.hpp"

#include "align/include/align_parameters.hpp"
#include "align/include/computeAlignments.hpp"
#include "align/include/parseCmdArgs.hpp"



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
     } else {
        robin_hood::unordered_flat_map< std::string, std::pair<skch::seqno_t, uint64_t> > seqName_to_seqCounterAndLen;
        skch::seqno_t seqCounter = 0;
        for(const auto &fileName : map_parameters.querySequences) {
            // check if there is a .fai
            std::string fai_name = fileName + ".fai";
            if (fs::exists(fai_name)) {
                // if so, process the .fai to determine our sequence length
                std::string line;
                std::ifstream in(fai_name.c_str());
                while (std::getline(in, line)) {
                    auto line_split = skch::CommonFunc::split(line, '\t');
                    const std::string seq_name = line_split[0];
                    const uint64_t seq_len = std::stoull(line_split[1]);
                    seqName_to_seqCounterAndLen[seq_name] = std::make_pair(seqCounter++,  seq_len);
                }
            } else {
                // if not, warn that this is expensive
                std::cerr << "[wfmash::align] WARNING, no .fai index found for " << fileName << ", reading the file to sort the mappings (slow)" << std::endl;
                for(const auto &fileName : map_parameters.querySequences)
                {
                    seqiter::for_each_seq_in_file(
						    fileName, {}, "", 
                            [&](const std::string& seq_name,
                                    const std::string& seq) {
                                seqName_to_seqCounterAndLen[seq_name] = std::make_pair(seqCounter++,  seq.length());
                            });
                }
            }
        }


        igzstream mappingListStream(map_parameters.outFileName.c_str());
        std::string mappingRecordLine;
        align::MappingBoundaryRow currentRecord;
        std::vector<align::MappingBoundaryRow> allReadMappings;

        while (!mappingListStream.eof()){
            std::getline(mappingListStream, mappingRecordLine);
            if( !mappingRecordLine.empty() ) {
                align::Aligner::parseMashmapRow(mappingRecordLine, currentRecord);

                allReadMappings.push_back(currentRecord);
            }
        }

        std::sort(allReadMappings.begin(), allReadMappings.end(), [&seqName_to_seqCounterAndLen](const align::MappingBoundaryRow &a, const align::MappingBoundaryRow &b)
        {
            return (seqName_to_seqCounterAndLen[a.qId].first < seqName_to_seqCounterAndLen[b.qId].first);
        });

        std::ofstream outstrm(align_parameters.mashmapPafFile);
        for(auto &e : allReadMappings)
        {
            outstrm << e.qId
            << "\t" << seqName_to_seqCounterAndLen[e.qId].second
            << "\t" << e.qStartPos
            << "\t" << e.qEndPos
            << "\t" << (e.strand == skch::strnd::FWD ? "+" : "-")
            << "\t" << e.refId
            << "\t" << seqName_to_seqCounterAndLen[e.refId].second
            << "\t" << e.rStartPos
            << "\t" << e.rEndPos
            << "\t" << 0
            << "\t" << std::max(e.rEndPos - e.rStartPos, e.qEndPos - e.qStartPos)
            << "\t" << 255
            << "\t" << "id:f:" << e.mashmap_estimated_identity
            << "\n";
        }
    }

    if (align_parameters.sam_format) {
        // Prepare SAM header
        std::ofstream outstrm(align_parameters.pafOutputFile);

        for(const auto &fileName : map_parameters.refSequences)
        {
            // check if there is a .fai
            std::string fai_name = fileName + ".fai";
            if (fs::exists(fai_name)) {
                // if so, process the .fai to determine our sequence length
                std::string line;
                std::ifstream in(fai_name.c_str());
                while (std::getline(in, line)) {
                    auto line_split = skch::CommonFunc::split(line, '\t');
                    const std::string seq_name = line_split[0];
                    const uint64_t seq_len = std::stoull(line_split[1]);
                    outstrm << "@SQ\tSN:" << seq_name << "\tLN:" << seq_len << "\n";
                }
            } else {
                // if not, warn that this is expensive
                std::cerr << "[wfmash::align] WARNING, no .fai index found for " << fileName << ", reading the file to prepare SAM header (slow)" << std::endl;
                seqiter::for_each_seq_in_file(
					    fileName, {}, "",
                        [&](const std::string& seq_name,
                                const std::string& seq) {
                            outstrm << "@SQ\tSN:" << seq_name << "\tLN:" << seq.length() << "\n";
                        });
            }





        }
        outstrm << "@PG\tID:wfmash\tPN:wfmash\tVN:0.1\tCL:wfmash\n";

        outstrm.close();
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
