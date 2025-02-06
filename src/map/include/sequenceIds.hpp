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

    seqno_t getSequenceId(const std::string& sequenceName) const {
        auto it = sequenceNameToId.find(sequenceName);
        if (it != sequenceNameToId.end()) {
            return it->second;
        }
        throw std::runtime_error("Sequence name not found: " + sequenceName);
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

        // Debug output to verify grouping
        for (size_t i = 0; i < metadata.size(); ++i) {
            std::cerr << "[DEBUG group assignment] seq " << metadata[i].name 
                      << " → group " << metadata[i].groupId << std::endl;
        }

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
        std::unordered_set<std::string> allowedQueryNames;
        std::unordered_set<std::string> allowedTargetNames;

        if (!queryList.empty()) readAllowedNames(queryList, allowedQueryNames);
        if (!targetList.empty()) readAllowedNames(targetList, allowedTargetNames);

        for (const auto& file : queryFiles) {
            readFAI(file, queryPrefixes, prefixDelim, allowedQueryNames, true);
        }
        for (const auto& file : targetFiles) {
            readFAI(file, targetPrefixes, prefixDelim, allowedTargetNames, false);
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
            return it->second;
        }
        seqno_t newId = metadata.size();
        sequenceNameToId[sequenceName] = newId;
        metadata.push_back(ContigInfo{sequenceName, length});
        return newId;
    }
};

} // namespace skch

#endif // SEQUENCE_ID_MANAGER_HPP
