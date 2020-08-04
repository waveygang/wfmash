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

#include "align/include/align_parameters.hpp"
#include "align/include/computeAlignments.hpp"

#include "yeet/include/parse_args.hpp"

//External includes
#include "common/args.hxx"

int main(int argc, char** argv) {
    /*
     * Make sure env variable MALLOC_ARENA_MAX is unset 
     * for efficient multi-thread execution
     */
    unsetenv((char *)"MALLOC_ARENA_MAX");

    args::ArgumentParser parser("edyeet: base-accurate alignments using edlib and mashmap2");
    args::HelpFlag help(parser, "help", "display this help menu", {'h', "help"});
    args::ValueFlag<uint64_t> thread_count(parser, "N", "use this many threads during parallel steps", {'t', "threads"});
    args::Positional<std::string> target_sequence_file(parser, "target", "alignment target or reference sequence file");
    args::ValueFlag<std::string> target_sequence_file_list(parser, "targets", "alignment target file list", {'L', "target-file-list"});
    args::PositionalList<std::string> query_sequences(parser, "queries", "query sequences");
    args::ValueFlag<std::string> query_sequence_file_list(parser, "queries", "alignment query file list", {'Q', "query-file-list"});
    // mashmap arguments
    args::Flag noSplit(parser, "no-split", "disable splitting of input sequences during mapping [enabled by default]", {'N',"no-split"});
    args::ValueFlag<float> map_pct_identity(parser, "%", "use this percent identity in the mashmap step [default: 85]", {'p', "map-pct-id"});
    args::ValueFlag<int> kmer_size(parser, "N", "kmer size <= 16 [default: 16]", {'k', "kmer"});
    args::ValueFlag<std::string> map_filter_mode(parser, "MODE", "filter mode for map step, either 'map', 'one-to-one', or 'none' [default: map]", {'f', "map-filter-mode"});
    args::ValueFlag<int> map_secondaries(parser, "N", "number of secondary mappings to retain in 'map' filter mode [default: 0]", {'n', "n-secondary"});
    // align parameters
    args::ValueFlag<std::string> align_input_paf(parser, "FILE", "derive precise alignments for this input PAF", {'i', "input-paf"});
    args::ValueFlag<float> align_pct_identity(parser, "%", "use this percent identity in the edlib step, if different than mashmap step [default: -p]", {'a', "align-pct-id"});
    args::ValueFlag<int> align_bandwidth(parser, "N", "maximum bandwidth for edlib alignment [default: 0 / computed from -p]", {'b', "align-bandwidth"});
    // 
    args::Flag keep_temp_files(parser, "", "keep intermediate files generated during mapping and alignment", {'T', "keep-temp"});
    args::Flag show_progress(parser, "show-progress", "write alignment progress to stderr", {'P', "show-progress"});
    args::Flag verbose_debug(parser, "verbose-debug", "enable verbose debugging", {'V', "verbose-debug"});

    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    if (argc==1) {
        std::cout << parser;
        return 1;
    }


    //CommandLineProcessing::ArgvParser cmd;

    // TODO:
    // if PAF input is given, skip mapping step and supply it directly to aligner
    // otherwise, if approximate mapping is requested, write map output directly to stdout
    // if base-level mapping is requested, write approximate mapping to file, and supply this to alignment

    // write map output into

    //Setup command line options
    //align::initCmdParser(cmd);
    //Parse command line arguements   

    //sketching and mapping parameters
    skch::Parameters map_parameters;
    align::Parameters align_parameters;
    yeet::parse_args(argc, argv, map_parameters, align_parameters);

    //parameters.refSequences.push_back(ref);

    //skch::parseandSave(argc, argv, cmd, parameters);

    auto t0 = skch::Time::now();

    //Build the sketch for reference
    skch::Sketch referSketch(map_parameters);

    std::chrono::duration<double> timeRefSketch = skch::Time::now() - t0;
    std::cerr << "INFO, skch::main, Time spent computing the reference index: " << timeRefSketch.count() << " sec" << std::endl;

    //Map the sequences in query file
    t0 = skch::Time::now();

    skch::Map mapper = skch::Map(map_parameters, referSketch);

    std::chrono::duration<double> timeMapQuery = skch::Time::now() - t0;
    std::cerr << "INFO, skch::main, Time spent mapping the query : " << timeMapQuery.count() << " sec" << std::endl;

    std::cerr << "INFO, skch::main, mapping results saved in : " << map_parameters.outFileName << std::endl;


    //Parse command line arguements   


    //align::parseandSave(argc, argv, cmd, parameters);   

    t0 = skch::Time::now();

    align::Aligner alignObj(align_parameters);

    std::chrono::duration<double> timeRefRead = skch::Time::now() - t0;
    std::cerr << "INFO, align::main, Time spent read the reference sequences: " << timeRefRead.count() << " sec" << std::endl;

    //Compute the alignments
    alignObj.compute();

    std::chrono::duration<double> timeAlign = skch::Time::now() - t0;
    std::cerr << "INFO, align::main, Time spent computing the aligment: " << timeAlign.count() << " sec" << std::endl;

    std::cerr << "INFO, align::main, alignment results saved in: " << align_parameters.samOutputFile << std::endl;
}
