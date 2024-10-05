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
    std::vector<std::string> idToSequenceName;
    std::vector<offset_t> idToSequenceLength;
    std::vector<std::string> querySequenceNames;
    std::vector<std::string> targetSequenceNames;

public:
    SequenceIdManager(const std::vector<std::string>& queryFiles,
                      const std::vector<std::string>& targetFiles,
                      const std::vector<std::string>& queryPrefixes,
                      const std::string& targetPrefix,
                      char prefixDelim,
                      const std::string& queryList = "",
                      const std::string& targetList = "") {
        populateFromFiles(queryFiles, targetFiles, queryPrefixes, targetPrefix, prefixDelim, queryList, targetList);
    }

    seqno_t getSequenceId(const std::string& sequenceName) const {
        auto it = sequenceNameToId.find(sequenceName);
        if (it != sequenceNameToId.end()) {
            return it->second;
        }
        throw std::runtime_error("Sequence name not found: " + sequenceName);
    }

    std::string getSequenceName(seqno_t id) const {
        if (id < static_cast<seqno_t>(idToSequenceName.size())) {
            return idToSequenceName[id];
        }
        throw std::runtime_error("Invalid sequence ID: " + std::to_string(id));
    }

    offset_t getSequenceLength(seqno_t id) const {
        if (id < static_cast<seqno_t>(idToSequenceLength.size())) {
            return idToSequenceLength[id];
        }
        throw std::runtime_error("Invalid sequence ID: " + std::to_string(id));
    }

    size_t size() const {
        return idToSequenceName.size();
    }

    const std::vector<std::string>& getQuerySequenceNames() const { return querySequenceNames; }
    const std::vector<std::string>& getTargetSequenceNames() const { return targetSequenceNames; }

private:
    void populateFromFiles(const std::vector<std::string>& queryFiles,
                           const std::vector<std::string>& targetFiles,
                           const std::vector<std::string>& queryPrefixes,
                           const std::string& targetPrefix,
                           char prefixDelim,
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
            readFAI(file, {targetPrefix}, prefixDelim, allowedTargetNames, false);
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
                 char prefixDelim,
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
        seqno_t newId = idToSequenceName.size();
        sequenceNameToId[sequenceName] = newId;
        idToSequenceName.push_back(sequenceName);
        idToSequenceLength.push_back(length);
        return newId;
    }
};

} // namespace skch

#endif // SEQUENCE_ID_MANAGER_HPP
