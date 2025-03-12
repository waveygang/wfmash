#ifndef SEQUENCE_ID_MANAGER_HPP
#define SEQUENCE_ID_MANAGER_HPP

#include <unordered_map>
#include <vector>
#include <string>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <unordered_set>
#include "base_types.hpp"

namespace skch {

class SequenceIdManager {
private:
    std::unordered_map<std::string, seqno_t> sequenceNameToId;
    std::vector<ContigInfo> metadata;
    std::vector<std::string> querySequenceNames;
    std::vector<std::string> targetSequenceNames;
    std::vector<std::string> allPrefixes;
    std::string prefixDelim;
    seqno_t nextId = 0;

public:
    SequenceIdManager(const std::vector<std::string>& queryFiles,
                      const std::vector<std::string>& targetFiles,
                      const std::vector<std::string>& queryPrefixes,
                      const std::vector<std::string>& targetPrefixes,
                      const std::string& prefixDelim,
                      const std::string& queryList = "",
                      const std::string& targetList = "")
        : prefixDelim(prefixDelim) {
        allPrefixes = queryPrefixes;
        allPrefixes.insert(allPrefixes.end(), targetPrefixes.begin(), targetPrefixes.end());
        populateFromFiles(queryFiles, targetFiles, queryPrefixes, targetPrefixes, prefixDelim, queryList, targetList);
        buildRefGroups();
    }
    
    // Export ID mapping information
    void exportIdMapping(std::ofstream& outStream) const {
        uint64_t mapSize = sequenceNameToId.size();
        outStream.write(reinterpret_cast<const char*>(&mapSize), sizeof(mapSize));
        
        for (const auto& [seqName, seqId] : sequenceNameToId) {
            uint64_t nameLength = seqName.size();
            outStream.write(reinterpret_cast<const char*>(&nameLength), sizeof(nameLength));
            outStream.write(seqName.c_str(), nameLength);
            outStream.write(reinterpret_cast<const char*>(&seqId), sizeof(seqId));
        }
        
        // Also export the next ID to use
        outStream.write(reinterpret_cast<const char*>(&nextId), sizeof(nextId));
    }
    
    // Import ID mapping information
    void importIdMapping(std::ifstream& inStream) {
        // Save original mappings in case we need to restore them
        auto originalMappings = sequenceNameToId;
        auto originalMetadata = metadata;
        auto originalNextId = nextId;
        
        try {
            // Clear current mappings to start fresh
            sequenceNameToId.clear();
            
            uint64_t mapSize = 0;
            inStream.read(reinterpret_cast<char*>(&mapSize), sizeof(mapSize));
            
            if (mapSize > 1000000) { // Sanity check - no realistic sequence set has millions of sequences
                throw std::runtime_error("Invalid mapping size in index file");
            }
            
            // First expand metadata to accommodate all IDs
            seqno_t maxId = 0;
            
            // Read all mappings from the index
            for (uint64_t i = 0; i < mapSize; ++i) {
                uint64_t nameLength = 0;
                inStream.read(reinterpret_cast<char*>(&nameLength), sizeof(nameLength));
                
                if (nameLength > 10000) { // Sanity check - no sequence name is this long
                    throw std::runtime_error("Invalid sequence name length in index file");
                }
                
                std::string seqName(nameLength, '\0');
                inStream.read(&seqName[0], nameLength);
                
                seqno_t seqId;
                inStream.read(reinterpret_cast<char*>(&seqId), sizeof(seqId));
                
                sequenceNameToId[seqName] = seqId;
                maxId = std::max(maxId, seqId);
            }
            
            // Resize metadata to accommodate all IDs
            metadata.resize(maxId + 1);
            
            // Update metadata for all mapped sequences
            for (const auto& [name, id] : sequenceNameToId) {
                if (id < metadata.size()) {
                    metadata[id].name = name;
                    
                    // Try to preserve sequence lengths from original metadata if available
                    auto origIt = originalMappings.find(name);
                    if (origIt != originalMappings.end() && origIt->second < originalMetadata.size()) {
                        metadata[id].len = originalMetadata[origIt->second].len;
                    }
                }
            }
            
            // Read the next ID
            seqno_t indexNextId;
            inStream.read(reinterpret_cast<char*>(&indexNextId), sizeof(indexNextId));
            
            // Ensure our nextId is at least as large as the one from the index
            nextId = std::max(indexNextId, maxId + 1);
            
        } catch (const std::exception& e) {
            // Restore original state on error
            sequenceNameToId = originalMappings;
            metadata = originalMetadata;
            nextId = originalNextId;
            std::cerr << "Error importing ID mappings: " << e.what() << std::endl;
            std::cerr << "Restored original ID mappings" << std::endl;
        }
    }

    seqno_t getSequenceId(const std::string& sequenceName) const {
        auto it = sequenceNameToId.find(sequenceName);
        if (it != sequenceNameToId.end()) {
            return it->second;
        }
        
        // Try a prefix-based lookup as a fallback
        for (const auto& [name, id] : sequenceNameToId) {
            if (name.find(sequenceName) == 0 || sequenceName.find(name) == 0) {
                std::cerr << "Warning: Using partial match for sequence '" << sequenceName 
                          << "' -> '" << name << "'" << std::endl;
                return id;
            }
        }
        
        // More detailed error message with available sequence names
        std::stringstream error_msg;
        error_msg << "Sequence name not found: '" << sequenceName << "'\nAvailable sequences:";
        int count = 0;
        for (const auto& [name, id] : sequenceNameToId) {
            if (count++ < 10) {
                error_msg << "\n  - '" << name << "' (ID: " << id << ")";
            }
        }
        if (count > 10) {
            error_msg << "\n  ... and " << (count - 10) << " more";
        }
        throw std::runtime_error(error_msg.str());
    }

    const ContigInfo& getContigInfo(seqno_t id) const {
        if (id < static_cast<seqno_t>(metadata.size())) {
            return metadata[id];
        }
        throw std::runtime_error("Invalid sequence ID: " + std::to_string(id));
    }

    const std::string& getSequenceName(seqno_t id) const {
        return getContigInfo(id).name;
    }

    const offset_t& getSequenceLength(seqno_t id) const {
        return getContigInfo(id).len;
    }

    size_t size() const {
        return metadata.size();
    }

    const std::vector<ContigInfo>& getMetadata() const {
        return metadata;
    }

    const std::vector<std::string>& getQuerySequenceNames() const { return querySequenceNames; }
    const std::vector<std::string>& getTargetSequenceNames() const { return targetSequenceNames; }

    int getRefGroup(seqno_t seqId) const {
        if (seqId < metadata.size()) {
            return metadata[seqId].groupId;
        }
        throw std::runtime_error("Invalid sequence ID: " + std::to_string(seqId));
    }

private:

    void buildRefGroups() {
        std::vector<std::tuple<std::string, size_t>> seqInfoWithIndex;
        size_t totalSeqs = metadata.size();

        for (size_t i = 0; i < totalSeqs; ++i) {
            seqInfoWithIndex.emplace_back(metadata[i].name, i);
        }

        std::sort(seqInfoWithIndex.begin(), seqInfoWithIndex.end());

        int currentGroup = 0;
        std::unordered_map<std::string, int> groupMap;

        for (const auto& [seqName, originalIndex] : seqInfoWithIndex) {
            std::string groupKey;

            if (!allPrefixes.empty()) {
                // Check if the sequence matches any of the specified prefixes
                // The original commented out version breaks in clang with an OpenMP capture error
                // auto it = std::find_if(allPrefixes.begin(), allPrefixes.end(), is_equal(prefix));
                    // [&seqName](const std::string& prefix) { return seqName.compare(0, prefix.length(), prefix) == 0; });

                for(auto prefix : allPrefixes) {
                    if (seqName.compare(0, prefix.length(), prefix) == 0) {
                        groupKey = prefix;
                        break;
                    }
                }
            }

            if (groupKey.empty() && !prefixDelim.empty()) {
                // Use prefix before last delimiter as group key
                size_t pos = seqName.rfind(prefixDelim);
                if (pos != std::string::npos) {
                    groupKey = seqName.substr(0, pos);
                }
            }

            if (groupKey.empty()) {
                // If no group key found, use the sequence name itself
                groupKey = seqName;
            }

            if (groupMap.find(groupKey) == groupMap.end()) {
                groupMap[groupKey] = ++currentGroup;
            }
            metadata[originalIndex].groupId = groupMap[groupKey];
        }
/*
        // Debug output to verify grouping
        for (size_t i = 0; i < metadata.size(); ++i) {
            std::cerr << "[DEBUG group assignment] seq " << metadata[i].name 
                      << " → group " << metadata[i].groupId << std::endl;
        }
*/
        if (totalSeqs == 0) {
            std::cerr << "[SequenceIdManager::buildRefGroups] ERROR: No sequences indexed!" << std::endl;
            exit(1);
        }
    }

    std::string getPrefix(const std::string& s) const {
        if (!prefixDelim.empty()) {
            size_t pos = s.find(prefixDelim);
            return (pos != std::string::npos) ? s.substr(0, pos) : s;
        }
        return s;
    }

    void populateFromFiles(const std::vector<std::string>& queryFiles,
                           const std::vector<std::string>& targetFiles,
                           const std::vector<std::string>& queryPrefixes,
                           const std::vector<std::string>& targetPrefixes,
                           const std::string& prefixDelim,
                           const std::string& queryList,
                           const std::string& targetList) {
        std::unordered_set<std::string> allowedTargetNames;
        std::unordered_set<std::string> allowedQueryNames;

        if (!targetList.empty()) readAllowedNames(targetList, allowedTargetNames);
        if (!queryList.empty()) readAllowedNames(queryList, allowedQueryNames);

        // Put target sequences first to keep their IDs stable with/without query sequences (useful for index generation)
        for (const auto& file : targetFiles) {
            readFAI(file, targetPrefixes, prefixDelim, allowedTargetNames, false);
        }
        for (const auto& file : queryFiles) {
            readFAI(file, queryPrefixes, prefixDelim, allowedQueryNames, true);
        }
    }

    void readAllowedNames(const std::string& listFile, std::unordered_set<std::string>& allowedNames) {
        std::ifstream file(listFile);
        std::string name;
        while (std::getline(file, name)) {
            allowedNames.insert(name);
        }
    }

    void readFAI(const std::string& fileName,
                 const std::vector<std::string>& prefixes,
                 const std::string& prefixDelim,
                 const std::unordered_set<std::string>& allowedNames,
                 bool isQuery) {
        std::string faiName = fileName + ".fai";
        std::ifstream faiFile(faiName);
        if (!faiFile.is_open()) {
            std::cerr << "Error: Unable to open FAI file: " << faiName << std::endl;
            exit(1);
        }

        std::string line;
        while (std::getline(faiFile, line)) {
            std::istringstream iss(line);
            std::string seqName;
            offset_t seqLength;
            iss >> seqName >> seqLength;

            bool prefixMatch = prefixes.empty() || std::any_of(prefixes.begin(), prefixes.end(),
                [&](const std::string& prefix) { return seqName.compare(0, prefix.size(), prefix) == 0; });

            if (prefixMatch && (allowedNames.empty() || allowedNames.find(seqName) != allowedNames.end())) {
                seqno_t seqId = addSequence(seqName, seqLength);
                if (isQuery) {
                    querySequenceNames.push_back(seqName);
                } else {
                    targetSequenceNames.push_back(seqName);
                }
            }
        }
    }

    seqno_t addSequence(const std::string& sequenceName, offset_t length) {
        auto it = sequenceNameToId.find(sequenceName);
        if (it != sequenceNameToId.end()) {
            // If we already have this sequence, update its length if needed
            if (metadata.size() <= it->second) {
                metadata.resize(it->second + 1);
                metadata[it->second] = ContigInfo{sequenceName, length};
            } else if (metadata[it->second].len != length) {
                metadata[it->second].len = length;
            }
            return it->second;
        }
        
        // Use nextId as the new sequence ID
        seqno_t newId = nextId++;
        sequenceNameToId[sequenceName] = newId;
        
        // Ensure metadata vector has enough capacity
        if (newId >= metadata.size()) {
            metadata.resize(newId + 1);
        }
        
        metadata[newId] = ContigInfo{sequenceName, length};
        return newId;
    }
};

} // namespace skch

#endif // SEQUENCE_ID_MANAGER_HPP
