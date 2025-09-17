/**
 * @file    externalSeeder.hpp
 * @brief   External seed input handler for PAF format
 * @author  Erik Garrison and Claude
 */

#ifndef EXTERNAL_SEEDER_HPP
#define EXTERNAL_SEEDER_HPP

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <optional>
#include <map>
#include <memory>
#include <atomic>

#include "map/include/base_types.hpp"
#include "map/include/map_parameters.hpp"
#include "map/include/sequenceIds.hpp"
#include "map/include/mappingFilter.hpp"
#include "map/include/mappingOutput.hpp"
#include "common/progress.hpp"

namespace skch {

// Helper structure to track PAF seed with query information
struct PAFSeed {
    MappingResult mapping;
    std::string queryName;
    offset_t queryLen;   // Store query sequence length from PAF
    offset_t queryEnd;   // Store query end since MappingResult only has start + blockLength
    offset_t refEnd;     // Store ref end for same reason
};

/**
 * @class     ExternalSeeder
 * @brief     Handler for processing external PAF seeds through wfmash filtering pipeline
 */
class ExternalSeeder {
public:
    /**
     * @brief Main entry point for processing external seeds
     * This injects PAF seeds directly into the wfmash filtering pipeline
     */
    static void processExternalSeeds(
        const Parameters& param,
        const std::string& seed_file,
        const SequenceIdManager& idManager,
        std::ostream& output_stream) {

        std::cerr << "[wfmash::externalSeeder] Reading external seeds from " << seed_file << std::endl;

        // Read and parse PAF seeds
        auto all_seeds = loadPAFSeeds(seed_file, idManager);

        if (all_seeds.empty()) {
            std::cerr << "[wfmash::externalSeeder] Warning: No valid seeds found in " << seed_file << std::endl;
            return;
        }

        std::cerr << "[wfmash::externalSeeder] Loaded " << all_seeds.size() << " seeds" << std::endl;

        // Group seeds by query sequence
        auto grouped = groupByQuery(all_seeds);

        std::cerr << "[wfmash::externalSeeder] Processing " << grouped.size() << " query sequences" << std::endl;

        // Process each query's seeds through the REAL wfmash filtering pipeline
        for (auto& [query_name, seeds] : grouped) {
            if (seeds.empty()) continue;

            // Convert PAFSeeds to MappingResults
            MappingResultsVector_t mappings;
            mappings.reserve(seeds.size());
            for (const auto& seed : seeds) {
                mappings.push_back(seed.mapping);
            }

            // Get query info
            seqno_t querySeqId = 0;
            offset_t queryLen = 0;

            // Get query length from the PAF seeds (they all have the same query)
            if (!seeds.empty()) {
                queryLen = seeds[0].queryLen;
            }

            // Try to get the sequence ID if it exists
            try {
                querySeqId = idManager.getSequenceId(query_name);
                // Also verify/update length from idManager if available
                offset_t idManagerLen = idManager.getSequenceLength(querySeqId);
                if (queryLen == 0) {
                    queryLen = idManagerLen;
                }
            } catch (...) {
                // Query not in ID manager, use what we have from PAF
            }

            // Create a progress meter for this query
            progress_meter::ProgressMeter progress(queryLen, "[wfmash::filter]", false);

            // Now inject these mappings into the ACTUAL filtering pipeline!
            // This is the key - we use the same filterSubsetMappings that the real mapper uses
            auto filteredResult = filterSubsetMappings(
                mappings,           // Our external seeds
                param,
                idManager,
                progress,
                querySeqId,
                queryLen
            );

            // Output the filtered mappings with proper chain info
            if (!filteredResult.merged_mappings.empty() || !filteredResult.unmerged_mappings.empty()) {
                // Use the merged mappings if merging was enabled, otherwise use unmerged
                auto& final_mappings = param.mergeMappings ?
                    filteredResult.merged_mappings : filteredResult.unmerged_mappings;
                auto& chain_info = filteredResult.chain_info;

                // Report the mappings using the standard output handler
                MappingOutput::reportReadMappings(
                    final_mappings,
                    chain_info,
                    query_name,
                    output_stream,
                    idManager,
                    param,
                    nullptr,  // no post-processing function
                    queryLen  // pass the query length
                );
            }
        }

        std::cerr << "[wfmash::externalSeeder] External seed processing complete" << std::endl;
    }

private:
    /**
     * @brief The actual filtering function that matches computeMap.hpp
     * This replicates the filterSubsetMappings function to work with external seeds
     */
    struct FilteredMappingsResult {
        MappingResultsVector_t merged_mappings;
        MappingResultsVector_t unmerged_mappings;
        ChainInfoVector_t chain_info;
    };

    static FilteredMappingsResult filterSubsetMappings(
        MappingResultsVector_t& mappings,
        const Parameters& param,
        const SequenceIdManager& idManager,
        progress_meter::ProgressMeter& progress,
        seqno_t querySeqId,
        offset_t queryLen) {

        FilteredMappingsResult result;

        if (mappings.empty()) return result;

        // Keep a copy of raw mappings for scaffold filtering
        MappingResultsVector_t rawMappings = mappings;

        // For external seeds, we skip the initial chaining (mergeMappingsInRange)
        // because FastGA seeds are already chained blocks
        // Instead, we just assign chain IDs to each mapping

        MappingsWithChains mappingsWithChains;
        if (param.mergeMappings) {
            // Even with mergeMappings enabled, for external seeds we don't merge
            // We just assign individual chain IDs
            mappingsWithChains.mappings = mappings;
            offset_t chain_id = 1;
            for (size_t i = 0; i < mappings.size(); i++) {
                ChainInfo ci;
                ci.chainId = chain_id++;
                ci.chainPos = 1;
                ci.chainLen = 1;
                mappingsWithChains.chainInfo.push_back(ci);
            }
        } else {
            // Same for non-merged path
            mappingsWithChains.mappings = mappings;
            offset_t chain_id = 1;
            for (size_t i = 0; i < mappings.size(); i++) {
                ChainInfo ci;
                ci.chainId = chain_id++;
                ci.chainPos = 1;
                ci.chainLen = 1;
                mappingsWithChains.chainInfo.push_back(ci);
            }
        }

        auto& workingMappings = mappingsWithChains.mappings;

        // Apply filtering based on parameters
        if (param.filterMode == filter::MAP || param.filterMode == filter::ONETOONE) {
            MappingResultsVector_t groupFilteredMappings;
            MappingFilterUtils::filterByGroup(
                workingMappings,
                groupFilteredMappings,
                param.numMappingsForSegment - 1,
                false,
                idManager,
                param,
                progress
            );
            workingMappings = std::move(groupFilteredMappings);
        }

        // Apply sparsification if needed
        MappingFilterUtils::sparsifyMappings(workingMappings, param);

        // Apply scaffold filtering if enabled and not in NONE filter mode
        if (param.scaffold_min_length > 0 && param.filterMode != filter::NONE) {
            // Create dummy progress trackers for scaffold filtering
            auto scaffold_progress = std::make_shared<progress_meter::ProgressMeter>(1, "[scaffold]", false);
            auto scaffold_total_work = std::make_shared<std::atomic<size_t>>(0);
            auto scaffold_completed_work = std::make_shared<std::atomic<size_t>>(0);

            MappingFilterUtils::filterByScaffolds(
                workingMappings,
                rawMappings,
                param,
                idManager,
                progress,
                querySeqId,
                queryLen,
                scaffold_progress,
                scaffold_total_work,
                scaffold_completed_work
            );
        }

        // Store results
        result.merged_mappings = workingMappings;
        result.unmerged_mappings = mappings;  // Original mappings
        result.chain_info = mappingsWithChains.chainInfo;

        return result;
    }

    /**
     * @brief Load and parse PAF seeds from file
     */
    static std::vector<PAFSeed> loadPAFSeeds(
        const std::string& seed_file,
        const SequenceIdManager& idManager) {

        std::vector<PAFSeed> all_seeds;
        std::ifstream infile(seed_file);

        if (!infile.is_open()) {
            if (seed_file == "-" || seed_file == "/dev/stdin") {
                // Try reading from stdin
                std::string line;
                while (std::getline(std::cin, line)) {
                    auto result = parsePAFLine(line, idManager);
                    if (result.has_value()) {
                        all_seeds.push_back(result.value());
                    }
                }
            } else {
                std::cerr << "[wfmash::externalSeeder] Error: Cannot open seed file " << seed_file << std::endl;
                exit(1);
            }
        } else {
            std::string line;
            size_t line_num = 0;
            while (std::getline(infile, line)) {
                line_num++;
                auto result = parsePAFLine(line, idManager);
                if (result.has_value()) {
                    all_seeds.push_back(result.value());
                } else {
                    std::cerr << "[wfmash::externalSeeder] Warning: Skipping invalid PAF line "
                              << line_num << std::endl;
                }
            }
            infile.close();
        }

        return all_seeds;
    }

    /**
     * @brief Parse a single PAF line into PAFSeed
     */
    static std::optional<PAFSeed> parsePAFLine(
        const std::string& line,
        const SequenceIdManager& idManager) {

        std::vector<std::string> fields;
        std::stringstream ss(line);
        std::string field;

        while (std::getline(ss, field, '\t')) {
            fields.push_back(field);
        }

        // PAF format requires at least 12 fields
        if (fields.size() < 12) {
            return std::nullopt;
        }

        PAFSeed result;
        MappingResult& mapping = result.mapping;

        // Parse required PAF fields
        result.queryName = fields[0];
        result.queryLen = std::stoull(fields[1]);  // Store query length
        mapping.queryStartPos = std::stoull(fields[2]);
        result.queryEnd = std::stoull(fields[3]);

        // Strand
        if (fields[4] == "+") {
            mapping.setStrand(strnd::FWD);
        } else if (fields[4] == "-") {
            mapping.setStrand(strnd::REV);
        } else {
            return std::nullopt; // Invalid strand
        }

        std::string target_name = fields[5];
        offset_t target_len = std::stoull(fields[6]);
        mapping.refStartPos = std::stoull(fields[7]);
        result.refEnd = std::stoull(fields[8]);

        // Get sequence ID from name
        try {
            mapping.refSeqId = idManager.getSequenceId(target_name);
        } catch (...) {
            std::cerr << "[wfmash::externalSeeder] Warning: Unknown target sequence '"
                      << target_name << "'" << std::endl;
            return std::nullopt;
        }

        // Calculate block length
        mapping.blockLength = result.refEnd - mapping.refStartPos;

        // Parse optional tags for identity/divergence
        mapping.setNucIdentity(0.9); // Default 90% identity

        for (size_t i = 12; i < fields.size(); i++) {
            if (fields[i].size() >= 5) {
                if (fields[i].substr(0, 5) == "dv:f:") {
                    // Divergence to identity conversion
                    float divergence = std::stof(fields[i].substr(5));
                    float identity = 1.0 - divergence;
                    mapping.setNucIdentity(identity);
                } else if (fields[i].substr(0, 5) == "id:f:") {
                    // Direct identity
                    float identity = std::stof(fields[i].substr(5));
                    mapping.setNucIdentity(identity);
                }
            }
        }

        // Set other required fields
        mapping.setKmerComplexity(1.0); // Full complexity for external seeds
        mapping.conservedSketches = 0; // Not applicable
        mapping.n_merged = 1;

        return result;
    }

    /**
     * @brief Group seeds by query sequence
     */
    static std::map<std::string, std::vector<PAFSeed>> groupByQuery(
        std::vector<PAFSeed>& all_seeds) {

        std::map<std::string, std::vector<PAFSeed>> grouped;

        for (auto& seed : all_seeds) {
            grouped[seed.queryName].push_back(seed);
        }

        return grouped;
    }
};

} // namespace skch

#endif // EXTERNAL_SEEDER_HPP