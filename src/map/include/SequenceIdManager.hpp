#ifndef SEQUENCE_ID_MANAGER_HPP
#define SEQUENCE_ID_MANAGER_HPP

#include <unordered_map>
#include <vector>
#include <string>
#include <stdexcept>
#include "base_types.hpp"

namespace skch {

class SequenceIdManager {
private:
    std::unordered_map<std::string, skch::seqno_t> sequenceNameToId;
    std::vector<std::string> idToSequenceName;

public:
    skch::seqno_t addSequence(const std::string& sequenceName) {
        auto it = sequenceNameToId.find(sequenceName);
        if (it != sequenceNameToId.end()) {
            return it->second;
        }
        skch::seqno_t newId = static_cast<skch::seqno_t>(idToSequenceName.size());
        sequenceNameToId[sequenceName] = newId;
        idToSequenceName.push_back(sequenceName);
        return newId;
    }

    skch::seqno_t getSequenceId(const std::string& sequenceName) const {
        auto it = sequenceNameToId.find(sequenceName);
        if (it != sequenceNameToId.end()) {
            return it->second;
        }
        throw std::runtime_error("Sequence name not found: " + sequenceName);
    }

    std::string getSequenceName(skch::seqno_t id) const {
        if (id < static_cast<skch::seqno_t>(idToSequenceName.size())) {
            return idToSequenceName[id];
        }
        throw std::runtime_error("Invalid sequence ID: " + std::to_string(id));
    }

    size_t size() const {
        return idToSequenceName.size();
    }
};

} // namespace skch

#endif // SEQUENCE_ID_MANAGER_HPP
