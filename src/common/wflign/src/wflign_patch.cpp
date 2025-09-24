#include <cstddef>
#include <chrono>
#include <cstdlib>
#include <string>
#include <atomic_image.hpp>
#include "rkmh.hpp"
#include "wflign_patch.hpp"
#include "wflign_git_version.hpp"

namespace wflign {

    std::tuple<std::string, uint64_t, uint64_t, uint64_t, uint64_t> trim_deletions(
        const std::string& cigar,
        uint64_t ref_start,
        uint64_t ref_end,
        uint64_t query_start,
        uint64_t query_end) {
        
        std::string trimmed_cigar;
        uint64_t new_ref_start = ref_start;
        uint64_t new_ref_end = ref_end;
        uint64_t new_query_start = query_start;
        uint64_t new_query_end = query_end;

        // Parse CIGAR string
        std::vector<std::pair<int, char>> cigar_ops;
        size_t i = 0;
        while (i < cigar.size()) {
            int count = 0;
            while (i < cigar.size() && std::isdigit(static_cast<unsigned char>(cigar[i]))) {
                count = count * 10 + (cigar[i] - '0');
                i++;
            }
            if (i < cigar.size()) {
                cigar_ops.push_back({count, cigar[i]});
                i++;
            }
        }

        // Process leading operations - keep insertions but remove deletions
        size_t leading_idx = 0;
        bool seen_non_id = false;
        std::vector<std::pair<int, char>> leading_insertions;
        
        // First pass: collect leading insertions and skip deletions
        while (leading_idx < cigar_ops.size()) {
            char op = cigar_ops[leading_idx].second;
            
            if (op == 'I') {
                // Save insertions to add later
                leading_insertions.push_back(cigar_ops[leading_idx]);
                leading_idx++;
            } else if (op == 'D') {
                // Skip deletions and adjust reference start
                new_ref_start += cigar_ops[leading_idx].first;
                leading_idx++;
            } else {
                // Stop at first non-I/D operation
                seen_non_id = true;
                break;
            }
        }
        
        // Process trailing operations - keep insertions but remove deletions
        size_t trailing_idx = cigar_ops.size() - 1;
        std::vector<std::pair<int, char>> trailing_insertions;
        
        // First pass from the end: collect trailing insertions and skip deletions
        // Only process if we have operations left after handling leading operations
        if (seen_non_id && leading_idx <= trailing_idx) {
            while (trailing_idx >= leading_idx) {
                char op = cigar_ops[trailing_idx].second;
                
                if (op == 'I') {
                    // Save insertions to add later
                    trailing_insertions.insert(trailing_insertions.begin(), cigar_ops[trailing_idx]);
                    if (trailing_idx == 0) break; // Prevent underflow
                    trailing_idx--;
                } else if (op == 'D') {
                    // Skip deletions and adjust reference end
                    new_ref_end -= cigar_ops[trailing_idx].first;
                    if (trailing_idx == 0) break; // Prevent underflow
                    trailing_idx--;
                } else {
                    // Stop at first non-I/D operation
                    break;
                }
            }
        }
        
        // Build new CIGAR string and count consumed bases
        // First add leading insertions
        for (const auto& op : leading_insertions) {
            trimmed_cigar += std::to_string(op.first) + op.second;
            // Update query position for insertions
            new_query_end += op.first;
        }
        
        // Add middle operations that weren't trimmed
        uint64_t ref_consumed = 0;
        uint64_t query_consumed = 0;
        
        for (size_t j = leading_idx; j <= trailing_idx && seen_non_id; j++) {
            trimmed_cigar += std::to_string(cigar_ops[j].first) + cigar_ops[j].second;
            
            // Count bases consumed by this operation
            switch(cigar_ops[j].second) {
                case 'M':
                case 'X':
                case '=':
                    ref_consumed += cigar_ops[j].first;
                    query_consumed += cigar_ops[j].first;
                    break;
                case 'D':
                case 'N':
                    ref_consumed += cigar_ops[j].first;
                    break;
                case 'I':
                    query_consumed += cigar_ops[j].first;
                    break;
                // S/H operations don't affect coordinates
            }
        }
        
        // Add trailing insertions
        for (const auto& op : trailing_insertions) {
            trimmed_cigar += std::to_string(op.first) + op.second;
            // Update query position for insertions
            query_consumed += op.first;
        }

        // Adjust coordinates based on consumed bases
        new_ref_end = new_ref_start + ref_consumed;
        new_query_end = new_query_start + query_consumed;

        return std::make_tuple(trimmed_cigar, new_ref_start, new_ref_end, new_query_start, new_query_end);
    }

    std::tuple<std::string, uint64_t, uint64_t, uint64_t, uint64_t> trim_indels(
        const std::string& cigar,
        uint64_t ref_start,
        uint64_t ref_end,
        uint64_t query_start,
        uint64_t query_end) {
    
        // Parse CIGAR string into operations
        std::vector<std::pair<int, char>> cigar_ops;
        size_t i = 0;
        while (i < cigar.size()) {
            int count = 0;
            while (i < cigar.size() && std::isdigit(static_cast<unsigned char>(cigar[i]))) {
                count = count * 10 + (cigar[i] - '0');
                i++;
            }
            if (i < cigar.size()) {
                cigar_ops.emplace_back(count, cigar[i]);
                i++;
            }
        }

        // Find first non-indel operation
        size_t start_idx = 0;
        uint64_t new_ref_start = ref_start;
        uint64_t new_query_start = query_start;

        while (start_idx < cigar_ops.size() && (cigar_ops[start_idx].second == 'I' || cigar_ops[start_idx].second == 'D')) {
            if (cigar_ops[start_idx].second == 'I') {
                new_query_start += cigar_ops[start_idx].first;
            } else { // 'D'
                new_ref_start += cigar_ops[start_idx].first;
            }
            start_idx++;
        }

        // Find last non-indel operation
        size_t end_idx = cigar_ops.size() - 1;
        uint64_t new_ref_end = ref_end;
        uint64_t new_query_end = query_end;
        
        if (start_idx < cigar_ops.size()) { // Only proceed if we have non-indel operations
            while (end_idx >= start_idx && (cigar_ops[end_idx].second == 'I' || cigar_ops[end_idx].second == 'D')) {
                if (cigar_ops[end_idx].second == 'I') {
                    new_query_end -= cigar_ops[end_idx].first;
                } else { // 'D'
                    new_ref_end -= cigar_ops[end_idx].first;
                }
                end_idx--;
            }
        }
        
        // Build trimmed CIGAR string
        std::string trimmed_cigar;
        uint64_t ref_consumed = 0;
        uint64_t query_consumed = 0;
        
        // Only build if we have operations to keep
        if (start_idx <= end_idx) {
            for (size_t j = start_idx; j <= end_idx; j++) {
                trimmed_cigar += std::to_string(cigar_ops[j].first) + cigar_ops[j].second;
                
                // Count bases consumed
                switch (cigar_ops[j].second) {
                    case 'M': case 'X': case '=':
                        ref_consumed += cigar_ops[j].first;
                        query_consumed += cigar_ops[j].first;
                        break;
                    case 'D': case 'N':
                        ref_consumed += cigar_ops[j].first;
                        break;
                    case 'I':
                        query_consumed += cigar_ops[j].first;
                        break;
                    // S/H operations don't affect coordinates
                }
            }
        }
    
        // Recalculate the end coordinates based on consumed bases
        new_ref_end = new_ref_start + ref_consumed;
        new_query_end = new_query_start + query_consumed;
    
        return std::make_tuple(trimmed_cigar, new_ref_start, new_ref_end, new_query_start, new_query_end);
    }

    // Process a compressed CIGAR string to get alignment metrics
    void process_compressed_cigar(
            const std::string& cigar_str,
            uint64_t& matches,
            uint64_t& mismatches, 
            uint64_t& insertions,
            uint64_t& inserted_bp,
            uint64_t& deletions,
            uint64_t& deleted_bp,
            uint64_t& refAlignedLength,
            uint64_t& qAlignedLength) {
        
        matches = mismatches = insertions = inserted_bp = deletions = deleted_bp = refAlignedLength = qAlignedLength = 0;
        
        
        // For a compressed CIGAR string like "50000=", we should get perfect identity
        if (cigar_str.find_first_not_of("0123456789") == cigar_str.length() - 1 && cigar_str.back() == '=') {
            matches = std::stoi(cigar_str.substr(0, cigar_str.length() - 1));
            refAlignedLength = matches;
            qAlignedLength = matches;
            return;
        }

        size_t pos = 0;
        while (pos < cigar_str.length()) {
            // Parse the length
            size_t next_pos = pos;
            while (next_pos < cigar_str.length() && isdigit(cigar_str[next_pos])) {
                next_pos++;
            }
            int length = std::stoi(cigar_str.substr(pos, next_pos - pos));
            char op = cigar_str[next_pos];
            
            switch (op) {
                case 'M':
                case '=':
                    matches += length;
                    refAlignedLength += length;
                    qAlignedLength += length;
                    break;
                case 'X':
                    mismatches += length;
                    refAlignedLength += length;
                    qAlignedLength += length;
                    break;
                case 'I':
                    insertions++;  // Count runs for gap-compressed
                    inserted_bp += length;  // Count individual bases for block
                    qAlignedLength += length;
                    break;
                case 'D':
                    deletions++;  // Count runs for gap-compressed
                    deleted_bp += length;  // Count individual bases for block
                    refAlignedLength += length;
                    break;
            }
            pos = next_pos + 1;
        }
    }

    void encodeOneStep(const char *filename, std::vector<unsigned char> &image, unsigned width, unsigned height) {
        //Encode the image
        unsigned error = lodepng::encode(filename, image, width, height);

        //if there's an error, display it
        if (error) std::cout << "encoder error " << error << ": " << lodepng_error_text(error) << std::endl;
    }

    namespace wavefront {

// accumulate alignment objects
// run the traceback determine which are part of the main chain
// order them and write them out
// needed--- 0-cost deduplication of alignment regions (how????)
//     --- trim the alignment back to the first 1/2 of the query
bool do_wfa_segment_alignment(
        const std::string& query_name,
        const char* query,
        std::vector<rkmh::hash_t>*& query_sketch,
        const uint64_t& query_length,
        const int64_t& j,
        const std::string& target_name,
        const char* target,
        std::vector<rkmh::hash_t>*& target_sketch,
        const uint64_t& target_length,
        const int64_t& i,
        const uint16_t& segment_length_q,
        const uint16_t& segment_length_t,
        const uint16_t& step_size,
        wflign_extend_data_t* extend_data,
        alignment_t& aln) {

    // if our i or j index plus segment length in the query or target is too long we'll make a memory access error and weird stuff will happen
    if (i + segment_length_t > target_length || j + segment_length_q > query_length) {
        // display function parameters
        std::cerr << "query_name: " << query_name << " query_length: " << query_length << " target_name: " << target_name << " target_length: " << target_length << std::endl;
        std::cerr << "i: " << i << " j: " << j << " segment_length_t: " << segment_length_t << " segment_length_q: " << segment_length_q << std::endl;
    }
    
    // first make the sketches if we haven't yet
    if (query_sketch == nullptr) {
        query_sketch = new std::vector<rkmh::hash_t>();
        *query_sketch = rkmh::hash_sequence(
                query + j, segment_length_q, extend_data->minhash_kmer_size, (uint64_t)((float)segment_length_q * extend_data->mash_sketch_rate));
        ++extend_data->num_sketches_allocated;
    }
    if (target_sketch == nullptr) {
        target_sketch = new std::vector<rkmh::hash_t>();        
        *target_sketch = rkmh::hash_sequence(
                target + i, segment_length_t, extend_data->minhash_kmer_size, (uint64_t)((float)segment_length_t * extend_data->mash_sketch_rate));
        ++extend_data->num_sketches_allocated;
    }

    // first check if our mash dist is inbounds
    const float mash_dist =
            rkmh::compare(*query_sketch, *target_sketch, extend_data->minhash_kmer_size);
    //std::cerr << "mash_dist is " << mash_dist << std::endl;

    // this threshold is set low enough that we tend to randomly sample wflambda
    // matrix cells for alignment the threshold is adaptive, based on the mash
    // distance of the mapping we are aligning we should obtain enough
    // alignments that we can still patch between them
    if (mash_dist > extend_data->max_mash_dist_to_evaluate) {
        // if it isn't, return false
        return false;
    } else {
        // if it is, we'll align
        const int max_score = (int)((float)std::max(segment_length_q, segment_length_t) * extend_data->inception_score_max_ratio);

        extend_data->wf_aligner->setMaxAlignmentSteps(max_score);
        const int status = extend_data->wf_aligner->alignEnd2End(
                target + i,segment_length_t,
                query + j,segment_length_q);

        aln.j = j;
        aln.i = i;

        aln.ok = (status == WF_STATUS_ALG_COMPLETED);

        // fill the alignment info if we aligned
        if (aln.ok) {
            aln.query_length = segment_length_q;
            aln.target_length = segment_length_t;

            /*
#ifdef VALIDATE_WFA_WFLIGN
            if (!validate_cigar(wf_aligner->cigar, query, target,
                        segment_length_q, segment_length_t, aln.j, aln.i)) {
                std::cerr << "cigar failure at alignment " << aln.j << " "
                          << aln.i << std::endl;
                unpack_display_cigar(wf_aligner->cigar, query,
                                     target, segment_length_q, segment_length_t,
                                     aln.j, aln.i);
                std::cerr << ">query" << std::endl
                          << std::string(query + j, segment_length_q)
                          << std::endl;
                std::cerr << ">target" << std::endl
                          << std::string(target + i, segment_length_t)
                          << std::endl;
                assert(false);
            }
#endif
             */

            wflign_edit_cigar_copy(*extend_data->wf_aligner,&aln.edit_cigar);

#ifdef VALIDATE_WFA_WFLIGN
            if (!validate_cigar(aln.edit_cigar, query, target, segment_length_q,
                        segment_length_t, aln.j, aln.i)) {
                std::cerr << "cigar failure after cigar copy in alignment "
                          << aln.j << " " << aln.i << std::endl;
                assert(false);
            }
#endif
        }

        return true;
    }
}

void do_wfa_patch_alignment(
        const char* query,
        const uint64_t& j,
        const uint64_t& query_length,
        const char* target,
        const uint64_t& i,
        const uint64_t& target_length,
        wfa::WFAlignerGapAffine2Pieces& wf_aligner,
        const wflign_penalties_t& convex_penalties,
        alignment_t& aln,
        alignment_t& rev_aln,
        const int64_t& chain_gap,
        const int& max_patching_score,
        const uint64_t& min_inversion_length) {

    const int max_score =
        max_patching_score ? max_patching_score :
            convex_penalties.gap_opening2 +
            (convex_penalties.gap_extension1 * std::min(
                    (int)chain_gap,
                    (int)std::max(target_length, query_length)
                )) + 64;

    /*
    std::cerr << "doing wfa patch alignment with parameters:"
              << " query_length " << query_length
              << " target_length " << target_length
              << " chain_gap " << chain_gap
              << " max_patching_score " << max_patching_score
              << " max_score " << max_score
              << std::endl
              << "query " << std::string(query + j, query_length)
              << " target " << std::string(target + i, target_length)
              << std::endl;
    */

    wf_aligner.setMaxAlignmentSteps(max_score);

    //int fwd_score = std::numeric_limits<int>::max();
    //int rev_score = std::numeric_limits<int>::max();
    
    const int status = wf_aligner.alignEnd2End(target + i, target_length, query + j, query_length);
    aln.ok = (status == WF_STATUS_ALG_COMPLETED);
    aln.is_rev = false;

    //std::cerr << "score is " << wf_aligner.getAlignmentScore() << std::endl;

    if (aln.ok) {
#ifdef VALIDATE_WFA_WFLIGN
        if (!validate_cigar(wf_aligner.cigar, query + j, target + i, query_length,
                            target_length, 0, 0)) {
            std::cerr << "cigar failure at alignment " << aln.j << " " << aln.i << std::endl;
            unpack_display_cigar(wf_aligner.cigar, query + j, target + i, query_length,
                                 target_length, 0, 0);
            std::cerr << ">query" << std::endl
                      << std::string(query + j, query_length) << std::endl;
            std::cerr << ">target" << std::endl
                      << std::string(target + i, target_length) << std::endl;
            assert(false);
        }
#endif

        wflign_edit_cigar_copy(wf_aligner, &aln.edit_cigar);
        aln.score = calculate_alignment_score(aln.edit_cigar, convex_penalties);
        //std::cerr << "forward score is " << fwd_score << std::endl;
    }

    if (query_length >= min_inversion_length && target_length >= min_inversion_length) {
        if (aln.ok) {
            wf_aligner.setMaxAlignmentSteps(std::ceil((double)aln.score * 0.9));
        }
        // Try reverse complement alignment
        std::string rev_comp_query = reverse_complement(std::string(query + j, query_length));
        const int rev_status = wf_aligner.alignEnd2End(target + i, target_length, rev_comp_query.c_str(), query_length);

        //auto rev_score = wf_aligner.getAlignmentScore();
        //rev_aln.ok = (rev_score > fwd_score && rev_status == WF_STATUS_ALG_COMPLETED);
        rev_aln.ok = (rev_status == WF_STATUS_ALG_COMPLETED);
        rev_aln.is_rev = true;

        if (rev_aln.ok) {
            wflign_edit_cigar_copy(wf_aligner, &rev_aln.edit_cigar);
            //std::cerr << "reverse complement alignment worked!" << std::endl;
#ifdef VALIDATE_WFA_WFLIGN
            if (!validate_cigar(rev_aln.edit_cigar, rev_comp_query.c_str(), target + i, query_length,
                                target_length, 0, 0)) {
                std::cerr << "cigar failure at reverse complement alignment " << j << " " << i << std::endl;
                unpack_display_cigar(rev_aln.edit_cigar, rev_comp_query.c_str(), target + i, query_length,
                                     target_length, 0, 0);
                std::cerr << ">query (reverse complement)" << std::endl
                          << rev_comp_query << std::endl;
                std::cerr << ">target" << std::endl
                          << std::string(target + i, target_length) << std::endl;
                assert(false);
            }
#endif

            rev_aln.score = calculate_alignment_score(rev_aln.edit_cigar, convex_penalties);
            //std::cerr << "reverse score is " << rev_score << std::endl;
            rev_aln.j = j;
            rev_aln.i = i;
            rev_aln.query_length = query_length;
            rev_aln.target_length = target_length;
        }
    } else {
        rev_aln.ok = false;
    }

    if (rev_aln.ok && rev_aln.score < aln.score) {
        rev_aln.ok = true;
        aln.ok = false;
#ifdef WFLIGN_DEBUG
        std::cerr << "got better score with reverse complement alignment" << std::endl
              << " query_length " << query_length
              << " target_length " << target_length
              << " chain_gap " << chain_gap
              << " max_patching_score " << max_patching_score
              << " max_score " << max_score
              << " fwd_score " << aln.score
              << " rev_score " << rev_aln.score
              << std::endl
              << "query " << std::string(query + j, query_length)
              << " target " << std::string(target + i, target_length)
              << std::endl;
#endif
    } else {
        rev_aln.ok = false;
    }

    aln.j = j;
    aln.i = i;
    aln.query_length = query_length;
    aln.target_length = target_length;
}

void erode_alignment(alignment_t& aln, int erode_k) {
    if (!aln.ok) return;

    std::vector<char> eroded_cigar;
    int match_count = 0;
    int non_match_count = 0;

    auto flush_matches = [&]() {
        if (match_count >= erode_k) {
            for (int i = 0; i < match_count; ++i) {
                eroded_cigar.push_back('M');
            }
        } else {
            for (int i = 0; i < match_count; ++i) {
                if (non_match_count > 0) {
                    eroded_cigar.push_back('D');
                    non_match_count--;
                } else {
                    eroded_cigar.push_back('I');
                }
            }
        }
        match_count = 0;
    };

    for (int i = aln.edit_cigar.begin_offset; i < aln.edit_cigar.end_offset; ++i) {
        char op = aln.edit_cigar.cigar_ops[i];
        if (op == 'M' || op == 'X') {
            if (non_match_count > 0) {
                flush_matches();
            }
            match_count++;
        } else {
            flush_matches();
            if (op == 'I') {
                eroded_cigar.push_back('I');
            } else if (op == 'D') {
                non_match_count++;
            }
        }
    }
    flush_matches();

    // Update the alignment
    free(aln.edit_cigar.cigar_ops);
    aln.edit_cigar.cigar_ops = (char*)malloc(eroded_cigar.size() * sizeof(char));
    memcpy(aln.edit_cigar.cigar_ops, eroded_cigar.data(), eroded_cigar.size() * sizeof(char));
    aln.edit_cigar.begin_offset = 0;
    aln.edit_cigar.end_offset = eroded_cigar.size();

    // Adjust query and target lengths
    int query_adjust = 0, target_adjust = 0;
    for (char op : eroded_cigar) {
        switch (op) {
            case 'M':
                query_adjust++;
                target_adjust++;
                break;
            case 'I':
                query_adjust++;
                break;
            case 'D':
                target_adjust++;
                break;
        }
    }
    aln.query_length = query_adjust;
    aln.target_length = target_adjust;
}

struct AlignmentBounds {
    int64_t query_start_offset;
    int64_t query_end_offset;
    int64_t target_start_offset;
    int64_t target_end_offset;
};

AlignmentBounds ok_find_alignment_bounds(const alignment_t& aln, const int& erode_k) {
    AlignmentBounds bounds;
    bounds.query_start_offset = 0;
    bounds.query_end_offset = 0;
    bounds.target_start_offset = 0;
    bounds.target_end_offset = 0;

    int64_t query_pos = 0;
    int64_t target_pos = 0;
    int match_count = 0;
    int reverse_match_count = 0;
    bool found_start = false;

    // Forward pass to find start bounds
    for (int i = aln.edit_cigar.begin_offset; i < aln.edit_cigar.end_offset; ++i) {
        char op = aln.edit_cigar.cigar_ops[i];
        switch (op) {
            case 'M':
            case '=':
                ++match_count;
                if (match_count >= erode_k && !found_start) {
                    bounds.query_start_offset = query_pos;
                    bounds.target_start_offset = target_pos;
                    found_start = true;
                }
                ++query_pos;
                ++target_pos;
                break;
            case 'X':
                match_count = 0;
                ++query_pos;
                ++target_pos;
                break;
            case 'I':
                match_count = 0;
                ++query_pos;
                break;
            case 'D':
                match_count = 0;
                ++target_pos;
                break;
        }
    }

    // Reverse pass to find end bounds
    query_pos = aln.query_length - 1;
    target_pos = aln.target_length - 1;
    for (int i = aln.edit_cigar.end_offset - 1; i >= aln.edit_cigar.begin_offset; --i) {
        char op = aln.edit_cigar.cigar_ops[i];
        switch (op) {
            case 'M':
            case '=':
                ++reverse_match_count;
                if (reverse_match_count >= erode_k) {
                    bounds.query_end_offset = query_pos + 1;
                    bounds.target_end_offset = target_pos + 1;
                    return bounds;
                }
                --query_pos;
                --target_pos;
                break;
            case 'X':
                reverse_match_count = 0;
                --query_pos;
                --target_pos;
                break;
            case 'I':
                reverse_match_count = 0;
                --query_pos;
                break;
            case 'D':
                reverse_match_count = 0;
                --target_pos;
                break;
        }
    }

    // If we didn't find end bounds, set them to the end of the alignment
    if (bounds.query_end_offset == 0) {
        bounds.query_end_offset = aln.query_length;
        bounds.target_end_offset = aln.target_length;
    }

    return bounds;
}

AlignmentBounds find_alignment_bounds(const alignment_t& aln, const int& erode_k) {
    AlignmentBounds bounds;
    bounds.query_start_offset = 0;
    bounds.query_end_offset = 0;
    bounds.target_start_offset = 0;
    bounds.target_end_offset = 0;

    int64_t query_pos = 0;
    int64_t target_pos = 0;
    int match_count = 0;
    bool found_start = false;
    bool found_end = false;

    // Forward pass
    for (int i = aln.edit_cigar.begin_offset; i < aln.edit_cigar.end_offset; ++i) {
        char op = aln.edit_cigar.cigar_ops[i];
        switch (op) {
            case 'M':
            case '=':
                ++match_count;
                if (match_count >= erode_k && !found_start) {
                    bounds.query_start_offset = query_pos;
                    bounds.target_start_offset = target_pos;
                    found_start = true;
                }
                ++query_pos;
                ++target_pos;
                break;
            case 'X':
                //match_count = 0;
                ++query_pos;
                ++target_pos;
                break;
            case 'I':
                //match_count = 0;
                ++query_pos;
                break;
            case 'D':
                //match_count = 0;
                ++target_pos;
                break;
        }
    }

    // Reverse pass
    query_pos = aln.query_length - 1;
    target_pos = aln.target_length - 1;
    match_count = 0;

    for (int i = aln.edit_cigar.end_offset - 1; i >= aln.edit_cigar.begin_offset; --i) {
        char op = aln.edit_cigar.cigar_ops[i];
        switch (op) {
            case 'M':
            case '=':
                ++match_count;
                if (match_count >= erode_k && !found_end) {
                    bounds.query_end_offset = query_pos + 1;
                    bounds.target_end_offset = target_pos + 1;
                    found_end = true;
                }
                --query_pos;
                --target_pos;
                break;
            case 'X':
                //match_count = 0;
                --query_pos;
                --target_pos;
                break;
            case 'I':
                //match_count = 0;
                --query_pos;
                break;
            case 'D':
                //match_count = 0;
                --target_pos;
                break;
        }
    }

    // If we didn't find bounds, set them to the full alignment length
    if (!found_start) {
        bounds.query_start_offset = 0;
        bounds.target_start_offset = 0;
    } else {
        // heuristic: subtract erode_k
        bounds.query_start_offset = std::max((int64_t)0, bounds.query_start_offset - erode_k);
        bounds.target_start_offset = std::max((int64_t)0, bounds.target_start_offset - erode_k);
    }
    if (!found_end) {
        bounds.query_end_offset = aln.query_length;
        bounds.target_end_offset = aln.target_length;
    } else {
        // heuristic: add erode_k
        bounds.query_end_offset = std::min((int64_t)aln.query_length, bounds.query_end_offset + erode_k);
        bounds.target_end_offset = std::min((int64_t)aln.target_length, bounds.target_end_offset + erode_k);
    }

    // Adjust bounds for reverse complement alignments
    if (aln.is_rev) {
        std::swap(bounds.query_start_offset, bounds.query_end_offset);
        bounds.query_start_offset = aln.query_length - bounds.query_start_offset;
        bounds.query_end_offset = aln.query_length - bounds.query_end_offset;
    }

    return bounds;
}

void trim_alignment(alignment_t& aln) {
    // Count leading indels
    int head_trim_q = 0, head_trim_t = 0;
    int head_trim_ops = 0;
    while (aln.edit_cigar.begin_offset + head_trim_ops < aln.edit_cigar.end_offset) {
        char op = aln.edit_cigar.cigar_ops[aln.edit_cigar.begin_offset + head_trim_ops];
        if (op != 'I' && op != 'D') break;
        if (op == 'I') head_trim_q++;
        if (op == 'D') head_trim_t++;
        head_trim_ops++;
    }

    // Count trailing indels
    int tail_trim_q = 0, tail_trim_t = 0;
    int tail_trim_ops = 0;
    while (aln.edit_cigar.end_offset - tail_trim_ops > aln.edit_cigar.begin_offset + head_trim_ops) {
        char op = aln.edit_cigar.cigar_ops[aln.edit_cigar.end_offset - 1 - tail_trim_ops];
        if (op != 'I' && op != 'D') break;
        if (op == 'I') tail_trim_q++;
        if (op == 'D') tail_trim_t++;
        tail_trim_ops++;
    }

    // If we found indels to trim, create new trimmed CIGAR
    if (head_trim_ops > 0 || tail_trim_ops > 0) {
        int new_len = aln.edit_cigar.end_offset - aln.edit_cigar.begin_offset - head_trim_ops - tail_trim_ops;
        char* new_cigar = (char*)malloc(new_len * sizeof(char));
        memcpy(new_cigar, 
               aln.edit_cigar.cigar_ops + aln.edit_cigar.begin_offset + head_trim_ops,
               new_len * sizeof(char));
        free(aln.edit_cigar.cigar_ops);
        aln.edit_cigar.cigar_ops = new_cigar;
        aln.edit_cigar.begin_offset = 0;
        aln.edit_cigar.end_offset = new_len;
    }

    // Adjust coordinates
    if (aln.is_rev) {
        aln.j += tail_trim_q;  // For reverse alignments, tail trim affects the start
    } else {
        aln.j += head_trim_q;
    }
    aln.i += head_trim_t;
    aln.query_length -= (head_trim_q + tail_trim_q);
    aln.target_length -= (head_trim_t + tail_trim_t);
}
        
std::vector<alignment_t> do_progressive_wfa_patch_alignment(
    const char* query,
    const uint64_t& query_start,
    const uint64_t& query_length,
    const char* target,
    const uint64_t& target_start,
    const uint64_t& target_length,
    wfa::WFAlignerGapAffine2Pieces& wf_aligner,
    const wflign_penalties_t& convex_penalties,
    const int64_t& chain_gap,
    const int& max_patching_score,
    const uint64_t& min_inversion_length,
    const int& erode_k) {

    std::vector<alignment_t> alignments;
    uint64_t current_query_start = query_start;
    uint64_t current_target_start = target_start;
    uint64_t remaining_query_length = query_length;
    uint64_t remaining_target_length = target_length;
    const uint64_t query_end = query_start + query_length;
    const uint64_t target_end = target_start + target_length;

    //std::cerr << "BEGIN do_progressive_wfa_patch_alignment: " << current_query_start << " " << remaining_query_length << " " << current_target_start << " " << remaining_target_length << std::endl;
    bool first = true;

    while (first || remaining_query_length >= min_inversion_length && remaining_target_length >= min_inversion_length) {
        first = false;
        //std::cerr << "do_progressive_wfa_patch_alignment: " << current_query_start << " " << remaining_query_length << " " << current_target_start << " " << remaining_target_length << std::endl;
        alignment_t aln, rev_aln;
        aln.is_rev = false;
        rev_aln.is_rev = true;
        do_wfa_patch_alignment(
            query,
            current_query_start,
            remaining_query_length,
            target,
            current_target_start,
            remaining_target_length,
            wf_aligner,
            convex_penalties,
            aln,
            rev_aln,
            chain_gap,
            max_patching_score,
            min_inversion_length);

        //std::cerr << "WFA fwd alignment: " << aln << std::endl;
        //std::cerr << "WFA rev alignment: " << rev_aln << std::endl;

        if (rev_aln.ok && (!aln.ok || rev_aln.score < aln.score)) {
            alignments.push_back(rev_aln);
        } else if (aln.ok) {
            alignments.push_back(aln);
            if (alignments.size() == 1) {
                break;
            }
        }

#ifdef VALIDATE_WFA_WFLIGN
        {
            std::string querystr = chosen_aln.is_rev ?
                reverse_complement(std::string(query + current_query_start, chosen_aln.query_length))
                : std::string(query + current_query_start, chosen_aln.query_length);
            std::string targetstr = std::string(target + current_target_start, chosen_aln.target_length);
            if (!validate_cigar(chosen_aln.edit_cigar, querystr.c_str(), targetstr.c_str(),
                                chosen_aln.query_length, chosen_aln.target_length, 0, 0)) {
                std::cerr << "before trim cigar failure at " << current_query_start << " " << current_target_start << std::endl;
                unpack_display_cigar(chosen_aln.edit_cigar, querystr.c_str(), targetstr.c_str(),
                                     chosen_aln.query_length, chosen_aln.target_length, 0, 0);
                std::cerr << ">query" << std::endl << querystr << std::endl
                          << ">target" << std::endl << targetstr << std::endl;
                assert(false);
                exit(1);
            }
        }
#endif

        if (alignments.empty()) {
            break;
        }
        auto& chosen_aln = alignments.back();
        auto bounds = find_alignment_bounds(chosen_aln, 7); // very light erosion of bounds on ends to avoid single-match starts and ends
        //std::cerr << "bounds: " << bounds.query_start_offset << " " << bounds.query_end_offset << " " << bounds.target_start_offset << " " << bounds.target_end_offset << std::endl;
                
#ifdef VALIDATE_WFA_WFLIGN
        {
            std::string querystr = chosen_aln.is_rev ?
                reverse_complement(std::string(query + current_query_start, chosen_aln.query_length))
                : std::string(query + current_query_start, chosen_aln.query_length);
            std::string targetstr = std::string(target + current_target_start, chosen_aln.target_length);
            if (!validate_cigar(chosen_aln.edit_cigar, querystr.c_str(), targetstr.c_str(),
                                chosen_aln.query_length, chosen_aln.target_length, 0, 0)) {
                std::cerr << "after trim cigar failure at " << current_query_start << " " << current_target_start << std::endl;
                unpack_display_cigar(chosen_aln.edit_cigar, querystr.c_str(), targetstr.c_str(),
                                     chosen_aln.query_length, chosen_aln.target_length, 0, 0);
                std::cerr << ">query" << std::endl << querystr << std::endl
                          << ">target" << std::endl << targetstr << std::endl;
                assert(false);
                exit(1);
            }
        }
#endif
        //alignments.push_back(chosen_aln);
        auto last_query_start = current_query_start;
        auto last_target_start = current_target_start;

        // Update the start positions and remaining lengths for the next iteration
        /*
        if (chosen_aln.is_rev) {
            current_query_start += bounds.query_start_offset;
        } else {
            current_query_start += bounds.query_end_offset;
        }
        current_target_start += bounds.target_end_offset;
        remaining_query_length = query_end - current_query_start;
        remaining_target_length = target_end - current_target_start;
        */
        // instead of the above, we want to progressively find the largest incomplete region to align
        // we should work with total area of the alignment over query and target to determine what to do
        // this will change the way we update tracking variables for the next iteration
        // basically: we take the larger of the two regions that will be at the start and end of the alignment
        // if there isn't any remaining target or query sequence on one end, obviously it's not an option
        // if there is, we should align the largest possible region that is still incomplete
        /// ....
        // Calculate the sizes of unaligned regions using only the bounds
        uint64_t left_query_size = bounds.query_start_offset;
        uint64_t right_query_size = remaining_query_length > bounds.query_end_offset ? remaining_query_length - bounds.query_end_offset : 0;
        uint64_t left_target_size = bounds.target_start_offset;
        uint64_t right_target_size = remaining_target_length > bounds.target_end_offset ? remaining_target_length - bounds.target_end_offset : 0;

        uint64_t max_left_size = std::max(left_query_size, left_target_size);
        uint64_t max_right_size = std::max(right_query_size, right_target_size);

        //std::cerr << "left_query_size: " << left_query_size << " left_target_size: " << left_target_size << std::endl
        //          << "right_query_size: " << right_query_size << " right_target_size: " << right_target_size << std::endl
        //          << "max_left_size: " << max_left_size << " max_right_size: " << max_right_size << std::endl;

        if (max_left_size >= max_right_size && max_left_size > 0) {
            // Align the left region
            remaining_query_length = left_query_size;
            remaining_target_length = left_target_size;
            // current_query_start and current_target_start remain unchanged
        } else if (max_right_size > 0) {
            // Align the right region
            current_query_start += bounds.query_end_offset;
            current_target_start += bounds.target_end_offset;
            remaining_query_length = right_query_size;
            remaining_target_length = right_target_size;
        } else {
            break; // No more regions to align
        }
    }
    //std::cerr << "END do_progressive_wfa_patch_alignment got " << alignments.size() << " alignments" << std::endl;
    return alignments;
}

void erode_head(std::vector<char>& unpatched, uint64_t& query_pos, uint64_t& target_pos, int erode_k) {
    int match_count = 0;
    auto it = unpatched.begin();
    uint64_t query_erased = 0, target_erased = 0;

    while (it != unpatched.end()) {
        if (*it == 'M' || *it == 'X') {
            match_count++;
            if (match_count >= erode_k) {
                break;
            }
            query_pos++;
            target_pos++;
            query_erased++;
            target_erased++;
        } else {
            //match_count = 0;
            if (*it == 'I') {
                query_pos++;
                query_erased++;
            }
            if (*it == 'D') {
                target_pos++;
                target_erased++;
            }
        }
        it++;
    }

    //std::cerr << "erode_head: eroded " << query_erased << " query and " << target_erased << " target" << std::endl;
    // Erase the eroded part
    unpatched.erase(unpatched.begin(), it);
}

void erode_tail(std::vector<char>& unpatched, uint64_t& query_end, uint64_t& target_end, int erode_k) {
    int match_count = 0;
    auto it = unpatched.rbegin();
    uint64_t q_offset = 0, t_offset = 0;

    while (it != unpatched.rend()) {
        if (*it == 'M' || *it == 'X') {
            match_count++;
            if (match_count >= erode_k) {
                break;
            }
            q_offset++;
            t_offset++;
        } else {
            //match_count = 0;
            if (*it == 'I') q_offset++;
            if (*it == 'D') t_offset++;
        }
        it++;
    }

    // Erase the eroded part
    //std::cerr << "erode_tail: eroded " << q_offset << " query and " << t_offset << " target" << std::endl;
    unpatched.erase(it.base(), unpatched.end());
    query_end -= q_offset;
    target_end -= t_offset;
}

void write_merged_alignment(
        std::ostream &out,
        const std::vector<alignment_t *> &trace,
        wfa::WFAlignerGapAffine2Pieces& wf_aligner,
        const wflign_penalties_t& convex_penalties,
        const bool& emit_md_tag,
        const bool& paf_format_else_sam,
        const bool& no_seq_in_sam,
        const char* query,
        const std::string& query_name,
        const uint64_t& query_total_length,
        const uint64_t& query_offset,
        const uint64_t& query_length,
        const bool& query_is_rev,
        const char* _target,
        const std::string& target_name,
        const uint64_t& target_total_length,
        const uint64_t& _target_offset,
        const uint64_t& _target_length,
        const float& min_identity,
        const uint64_t& min_alignment_length,
        const float& min_block_identity,
#ifdef WFA_PNG_TSV_TIMING
        const long& elapsed_time_wflambda_ms,
        const uint64_t& num_alignments,
        const uint64_t& num_alignments_performed,
#endif
        const float& mashmap_estimated_identity,
        const uint64_t& wflign_max_len_major,
        const uint64_t& wflign_max_len_minor,
        const int& erode_k,
        const int64_t& chain_gap,
        const int& max_patching_score,
        const uint64_t& min_inversion_length,
        const int& min_wf_length,
        const int& max_dist_threshold,
#ifdef WFA_PNG_TSV_TIMING
        const std::string* prefix_wavefront_plot_in_png,
        const uint64_t& wfplot_max_size,
        const bool& emit_patching_tsv,
        std::ostream* out_patching_tsv,
#endif
        const bool& with_endline) {

    int64_t target_pointer_shift = 0;

    // we need to get the start position in the query and target
    // then run through the whole alignment building up the cigar
    // finally emitting it
    // our final cigar
    //
    // std::string cigarstr;
    uint64_t matches = 0;
    uint64_t mismatches = 0;
    uint64_t insertions = 0;
    uint64_t inserted_bp = 0;
    uint64_t deletions = 0;
    uint64_t deleted_bp = 0;
    uint64_t query_start = 0;
    uint64_t target_start = 0;
    uint64_t total_query_aligned_length = 0;
    uint64_t total_target_aligned_length = 0;
    uint64_t query_end = 0;
    uint64_t target_end = 0;
    //uint64_t total_score = 0;

    char* target = (char*)_target;
    uint64_t target_offset = _target_offset;
    uint64_t target_length = _target_length;

    // double mash_dist_sum = 0;
    uint64_t ok_alns = 0;

    auto start_time = std::chrono::steady_clock::now();

    //std::cerr << "target_offset: " << target_offset << " query_offset: " << query_offset << std::endl;
    //std::cerr << "target_length: " << target_length << " query_length: " << query_length << std::endl;
    //std::cerr << "target_total_length: " << target_total_length << " query_total_length: " << query_total_length << std::endl;

#ifdef WFLIGN_DEBUG
    std::cerr << "[wflign::wflign_affine_wavefront] processing traceback"
      << std::endl;
#endif
    // write trace into single cigar vector
    std::vector<char> tracev;
    std::vector<alignment_t> multi_patch_alns;
    {
        // patch: walk the cigar, patching directly when we have simultaneous
        // gaps in query and ref and adding our results to the final trace as we
        // go

#define MAX_NUM_INDELS_TO_LOOK_AT 2
        auto distance_close_big_enough_indels =	
                [](const uint32_t indel_len, auto iterator,	
                   const std::vector<char> &trace,
                   const uint16_t&max_dist_to_look_at) {	
                    const uint32_t min_indel_len_to_find = indel_len / 3;	

                    // std::cerr << "min_indel_len_to_find " <<	
                    // min_indel_len_to_find << std::endl; std::cerr <<	
                    // "max_dist_to_look_at " << max_dist_to_look_at << std::endl;	

                    auto q = iterator;	

                    uint8_t num_indels_to_find = MAX_NUM_INDELS_TO_LOOK_AT;	
                    uint32_t curr_size_close_indel = 0;	
                    int32_t dist_close_indels = 0;	

                    // std::cerr << "indel_len " << indel_len << std::endl;	
                    // std::cerr << "min_indel_len_to_find " << min_indel_len_to_find << std::endl;	
                    // std::cerr << "max_dist_to_look_at " << max_dist_to_look_at << std::endl;	

                    while (q != trace.end() &&	
                           dist_close_indels < max_dist_to_look_at) {	
                        curr_size_close_indel = 0;	
                        while (q != trace.end() && (*q == 'I' || *q == 'D')) {	
                            ++curr_size_close_indel;	

                            ++dist_close_indels;	
                            ++q;	
                        }	
                        // std::cerr << "\t\tcurr_size_close_indel " <<	
                        // curr_size_close_indel << std::endl;	
                        if (curr_size_close_indel >= min_indel_len_to_find) {	
                            // std::cerr << "\t\tnum_indels_to_find " <<	
                            // (uint16_t)num_indels_to_find << std::endl;	
                            if (--num_indels_to_find == 0) {	
                                break;	
                            }	
                        }	

                        while (q != trace.end() &&	
                               (dist_close_indels < max_dist_to_look_at) &&	
                               *q != 'I' && *q != 'D') {	
                            ++dist_close_indels;	
                            ++q;	
                        }	
                    }	

                    // std::cerr << "dist_close_indels " << dist_close_indels << std::endl;	
                    // std::cerr << "num_indels_found " << MAX_NUM_INDELS_TO_LOOK_AT - num_indels_to_find << std::endl;	

                    return num_indels_to_find < MAX_NUM_INDELS_TO_LOOK_AT	
                           ? dist_close_indels	
                           : -1;	
                };

        auto patching = [&query, &query_name, &query_length, &query_start,
                         &query_offset, &query_end, &target, &target_name,
                         &target_length, &target_start, &target_offset,
                         &target_total_length, &target_end,
                         &target_pointer_shift,
                         &wflign_max_len_major,
                         &wflign_max_len_minor,
                         &distance_close_big_enough_indels, &min_wf_length,
                         &max_dist_threshold, &wf_aligner,
                         &multi_patch_alns,
                         &convex_penalties,
                         &chain_gap, &max_patching_score, &min_inversion_length, &erode_k,
                         &query_total_length  // Add this line to capture query_total_length
#ifdef WFA_PNG_TSV_TIMING
                         ,&emit_patching_tsv,
                         &out_patching_tsv
#endif
            ](std::vector<char> &unpatched,
              std::vector<char> &patched,
              const uint16_t &min_wfa_head_tail_patch_length,
              const uint16_t &min_wfa_patch_length,
              const uint16_t &max_dist_to_look_at,
              bool save_multi_patch_alns) {

            auto q = unpatched.begin();

            uint64_t query_delta = 0;
            uint64_t target_delta = 0;

            bool got_alignment = false;

            // trim back small matches at the start
            erode_head(unpatched, query_start, target_start, 7);

            uint64_t query_pos = query_start;
            uint64_t target_pos = target_start;

            // Trim spurious matches off the tail of the alignment
            erode_tail(unpatched, query_end, target_end, 7);

            // Head patching
            if (query_start > 0 || target_start > 0) {
                // Calculate how far we need to shift to cover the query, and how far we can safely shift
                int64_t needed_shift = std::max(0L, static_cast<int64_t>(query_start) - static_cast<int64_t>(target_start));
                int64_t max_safe_shift = std::min(
                    static_cast<int64_t>(wflign_max_len_minor),
                    static_cast<int64_t>(target_offset)
                    );
    
                // Take the minimum of what we need and what's safe
                int64_t actual_shift = std::min(needed_shift, max_safe_shift);

                // Adjust the target pointer, offset, and length
                target -= actual_shift;
                target_offset -= actual_shift;
                target_length += actual_shift;
                target_end += actual_shift; // hmm
                target_pos += actual_shift;
                target_start += actual_shift;

                alignment_t head_aln;
                alignment_t head_rev_aln;
                do_wfa_patch_alignment(
                    query,
                    0,
                    query_start,
                    target,
                    0,  // Start from the beginning of the adjusted target
                    target_start,
                    wf_aligner,
                    convex_penalties,
                    head_aln,
                    head_rev_aln,
                    chain_gap,
                    max_patching_score,
                    min_inversion_length);
                
                if (head_aln.ok) {
                    //std::cerr << "head_aln: " << head_aln.score << std::endl;
                    // Prepend the head alignment to the main alignment
                    for (int i = head_aln.edit_cigar.begin_offset; i < head_aln.edit_cigar.end_offset; ++i) {
                        //std::cerr << head_aln.edit_cigar.cigar_ops[i];
                        patched.push_back(head_aln.edit_cigar.cigar_ops[i]);
                    }
                    //std::cerr << std::endl;
                } else {
                    // push back I and D to fill the gap
                    for (uint64_t i = 0; i < query_start; ++i) {
                        patched.push_back('I');
                    }
                    for (uint64_t i = 0; i < target_start; ++i) {
                        patched.push_back('D');
                    }
                }
                query_start = 0;
                target_start = 0;
            }

            // Patching in the middle
            while (q != unpatched.end()) {
                // get to the first matchn
                while (q != unpatched.end() && (*q == 'M' || *q == 'X')) {
                    /*
                std::cerr << "q: " << query[query_pos] << " "
                          << "t: " << target[target_pos - target_pointer_shift]
                          << std::endl;
                */
                    if (query_pos >= query_length ||
                        target_pos >= target_length) {
                        std::cerr << "[wflign::wflign_affine_wavefront] "
                                     "corrupted traceback (out of bounds) for "
                                  << query_name << " " << query_offset << " "
                                  << target_name << " " << target_offset
                                  << std::endl;
                        exit(1);
                    }

                    if (*q == 'M') {
                        if (query[query_pos] !=
                            target[target_pos - target_pointer_shift]) {
                            std::cerr << "[wflign::wflign_affine_wavefront] "
                                         "corrupted traceback (M, but there is "
                                         "a mismatch) for "
                                      << query_name << " " << query_offset
                                      << " " << target_name << " "
                                      << target_offset << std::endl;
                            exit(1);
                        }
                    } else {
                        if (query[query_pos] ==
                            target[target_pos - target_pointer_shift]) {
                            std::cerr << "[wflign::wflign_affine_wavefront] "
                                         "corrupted traceback (X, but there is "
                                         "a match) for "
                                      << query_name << " " << query_offset
                                      << " " << target_name << " "
                                      << target_offset << std::endl;
                            exit(1);
                        }
                    }

                    patched.push_back(*q);
                    ++query_pos;
                    ++target_pos;
                    ++q;
                }

                // how long a gap?
                while (q != unpatched.end() && *q == 'I') {
                    ++query_delta;
                    ++q;
                }
                while (q != unpatched.end() && *q == 'D') {
                    ++target_delta;
                    ++q;
                }

                // how long was our last gap?
                // if it's long enough, patch it
                int32_t size_region_to_repatch = 0;
                // unused!

                {
                    got_alignment = false;

                    if ((size_region_to_repatch > 0 ||
                         (query_delta > 0 && target_delta > 0) ||
                         (query_delta > 2 || target_delta > 2) &&
                                 (query_delta < wflign_max_len_major &&
                                  target_delta < wflign_max_len_major) &&
                                 (query_delta < wflign_max_len_minor ||
                                  target_delta < wflign_max_len_minor))) {

                        int32_t distance_close_indels = 
                            (query_delta > 10 || target_delta > 10) ?	
                            distance_close_big_enough_indels(std::max(query_delta, target_delta), q, unpatched, max_dist_to_look_at) :	
                            -1;

                        // Trigger the patching if there is a dropout
                        // (consecutive Is and Ds) or if there is a close and
                        // big enough indel forward
                        if (size_region_to_repatch > 0 ||
                            (query_delta > 0 && target_delta > 0) ||
                            (query_delta > 2 || target_delta > 2) ||
                            distance_close_indels > 0) {
#ifdef WFLIGN_DEBUG
                            // std::cerr << "query_delta " << query_delta <<
                            // "\n"; std::cerr << "target_delta " << target_delta
                            // << "\n"; std::cerr << "distance_close_indel " <<
                            // distance_close_indel << "\n";

                            std::cerr << "[wflign::wflign_affine_wavefront] "
                                         "patching in "
                                      << query_name << " " << query_offset
                                      << " @ " << query_pos << " - "
                                      << query_delta << " " << target_name
                                      << " " << target_offset << " @ "
                                      << target_pos << " - " << target_delta
                                      << std::endl;
#endif

                            // if we are continuing a patch, we can't nibble
                            // backward too much to avoid the risk of going in
                            // endless loop
                            if (size_region_to_repatch > 0) {
                                // nibble backward
                                while (!patched.empty() &&
                                       size_region_to_repatch > 0) {
                                    const auto &c = patched.back();
                                    switch (c) {
                                        case 'M':
                                        case 'X':
                                            --query_pos;
                                            --target_pos;
                                            ++query_delta;
                                            ++target_delta;
                                            break;
                                        case 'I':
                                            ++query_delta;
                                            --query_pos;
                                            break;
                                        case 'D':
                                            ++target_delta;
                                            --target_pos;
                                            break;
                                        default:
                                            break;
                                    }
                                    patched.pop_back();
                                    --size_region_to_repatch;
                                }
                            } else {
                                // nibble backward if we're below the correct
                                // length
                                while (
                                        !patched.empty() &&
                                        (query_delta < (min_wfa_patch_length / 2) ||
                                         target_delta <
                                                 (min_wfa_patch_length / 2))) {
                                    const auto &c = patched.back();
                                    switch (c) {
                                        case 'M':
                                        case 'X':
                                            --query_pos;
                                            --target_pos;
                                            ++query_delta;
                                            ++target_delta;
                                            break;
                                        case 'I':
                                            ++query_delta;
                                            --query_pos;
                                            break;
                                        case 'D':
                                            ++target_delta;
                                            --target_pos;
                                            break;
                                        default:
                                            break;
                                    }
                                    patched.pop_back();
                                }
                            }

                            // nibble forward if we're below the correct length
                            while (q != unpatched.end() &&
                                   (query_delta < min_wfa_patch_length ||
                                    target_delta < min_wfa_patch_length)) {
                                const auto &c = *q++;
                                switch (c) {
                                    case 'M':
                                    case 'X':
                                        ++query_delta;
                                        ++target_delta;
                                        break;
                                    case 'I':
                                        ++query_delta;
                                        break;
                                    case 'D':
                                        ++target_delta;
                                        break;
                                    default:
                                        break;
                                }

                                --distance_close_indels;
                            }

                            // Nibble until the close, big enough indel is	
                            // reached Important when the patching can't be	
                            // computed correctly without including the next	
                            // indel	
                            while (q != unpatched.end() &&	
                                   distance_close_indels > 0) {	
                                const auto &c = *q++;	
                                switch (c) {	
                                    case 'M':	
                                    case 'X':	
                                        ++query_delta;	
                                        ++target_delta;	
                                        break;	
                                    case 'I':	
                                        ++query_delta;	
                                        break;	
                                    case 'D':	
                                        ++target_delta;	
                                        break;	
                                    default:	
                                        break;	
                                }	

                                --distance_close_indels;	
                            }

                            // check forward if there are other Is/Ds to merge
                            // in the current patch
                            while (q != unpatched.end() &&
                                   (*q == 'I' || *q == 'D') &&
                                   ((query_delta < wflign_max_len_major &&
                                     target_delta < wflign_max_len_major) &&
                                    (query_delta < wflign_max_len_minor ||
                                     target_delta < wflign_max_len_minor))) {
                                const auto &c = *q++;
                                if (c == 'I') {
                                    ++query_delta;
                                } else {
                                    ++target_delta;
                                }
                            }

                            // check backward if there are other Is/Ds to merge
                            // in the current patch it will eventually nibble
                            // the Is/Ds left from the last patch
                            while (!patched.empty() &&
                                   (patched.back() == 'I' ||
                                    patched.back() == 'D') &&
                                   ((query_delta < wflign_max_len_major &&
                                     target_delta < wflign_max_len_major) &&
                                    (query_delta < wflign_max_len_minor ||
                                     target_delta < wflign_max_len_minor))) {
                                const auto &c = patched.back();
                                if (c == 'I') {
                                    ++query_delta;
                                    --query_pos;
                                } else {
                                    ++target_delta;
                                    --target_pos;
                                }
                                patched.pop_back();
                            }

                            size_region_to_repatch = 0;
                            {
                                // WFA is only global
                                auto patch_alignments = do_progressive_wfa_patch_alignment(
                                        query,
                                        query_pos,
                                        query_delta,
                                        target - target_pointer_shift,
                                        target_pos,
                                        target_delta,
                                        wf_aligner,
                                        convex_penalties,
                                        chain_gap,
                                        max_patching_score,
                                        min_inversion_length,
                                        erode_k);
                                if (patch_alignments.size() == 1
                                    && patch_alignments.front().ok
                                    && !patch_alignments.front().is_rev) {
                                    got_alignment = true;
                                    auto& patch_aln = patch_alignments.front();
                                    const int start_idx =
                                        patch_aln.edit_cigar.begin_offset;
                                    const int end_idx =
                                        patch_aln.edit_cigar.end_offset;
                                    for (int i = start_idx; i < end_idx; i++) {
                                        patched.push_back(patch_aln.edit_cigar.cigar_ops[i]);
                                    }
                                } else if (save_multi_patch_alns) {
                                    for (auto& aln : patch_alignments) {
                                        trim_alignment(aln);
                                        multi_patch_alns.push_back(aln);
                                    }
                                }

#ifdef WFA_PNG_TSV_TIMING
                                if (emit_patching_tsv) {
                                    for (auto& aln : patch_alignments) {
                                        *out_patching_tsv
                                            << query_name << "\t" << query_pos << "\t" << query_pos + query_delta << "\t"
                                            << target_name << "\t" << (target_pos - target_pointer_shift) << "\t" << (target_pos - target_pointer_shift + target_delta) << "\t"
                                            << aln.ok << std::endl;
                                    }
                                }
#endif
                            }
                        }
                    }

                    // add in stuff if we didn't align
                    if (!got_alignment) {
                        for (uint64_t i = 0; i < query_delta; ++i) {
                            patched.push_back('I');
                        }
                        for (uint64_t i = 0; i < target_delta; ++i) {
                            patched.push_back('D');
                        }
                    }

                    // std::cerr << "query_delta " << query_delta << std::endl;
                    // std::cerr << "target_delta " << target_delta <<
                    // std::endl;
                    query_pos += query_delta;
                    target_pos += target_delta;

                    query_delta = 0;
                    target_delta = 0;
                }
            }

            // trim Is and Ds from the end of the alignment trace, adjusting the query_pos and target_pos as needed
            while (!patched.empty() && (patched.back() == 'I' || patched.back() == 'D')) {
                if (patched.back() == 'I') {
                    --query_pos;
                } else {
                    --target_pos;
                }
                patched.pop_back();
            }

            // Tail patching
            if (query_pos < query_length || target_pos < target_length) {
                // Calculate how much additional target sequence we need
                int64_t needed_extension = std::max(0L, static_cast<int64_t>(query_length - query_pos) - static_cast<int64_t>(target_length - target_pos));
    
                // Calculate how much we can safely extend
                int64_t max_safe_extension = std::min(
                    static_cast<int64_t>(wflign_max_len_minor),
                    static_cast<int64_t>(target_total_length - (target_offset + target_length))
                    );

                // Take the minimum of what we need and what's safe
                int64_t actual_extension = std::min(needed_extension, max_safe_extension);

                alignment_t tail_aln;
                alignment_t tail_rev_aln;
                do_wfa_patch_alignment(
                    query,
                    query_pos,
                    query_length - query_pos,
                    target,
                    target_pos,
                    (target_length - target_pos) + actual_extension,
                    wf_aligner,
                    convex_penalties,
                    tail_aln,
                    tail_rev_aln,
                    chain_gap,
                    max_patching_score,
                    min_inversion_length);

                if (tail_aln.ok) {
                    // Append the tail alignment to the main alignment
                    for (int i = tail_aln.edit_cigar.begin_offset; i < tail_aln.edit_cigar.end_offset; ++i) {
                        patched.push_back(tail_aln.edit_cigar.cigar_ops[i]);
                    }
                    query_pos = query_length;
                    target_pos = target_length;

                    // Calculate new ends relative to the segment being aligned
                    uint64_t new_query_end = query_length;
                    uint64_t new_target_end = target_length + actual_extension;

                    if (query_offset + new_query_end > query_total_length || target_offset + new_target_end > target_total_length) {
                        std::cerr << "[wfmash::patch] Warning: Alignment extends beyond sequence bounds. Truncating." << std::endl;
                    }

                    // Adjust query_end and target_end, ensuring we don't exceed the segment lengths
                    query_end = std::min(new_query_end, query_length);
                    target_end = std::min(new_target_end, target_length + actual_extension);
                }
            }

#ifdef WFLIGN_DEBUG
            std::cerr << "[wflign::wflign_affine_wavefront] got unsorted "
                 "patched traceback: ";
            for (auto c : patched) {
                std::cerr << c;
            }
            std::cerr << std::endl;
#endif
        };

        {
            std::vector<char> erodev;
            {
                std::vector<char> rawv;

                // copy
#ifdef WFLIGN_DEBUG
                std::cerr
                    << "[wflign::wflign_affine_wavefront] copying traceback"
                    << std::endl;
#endif
                for (auto x = trace.rbegin(); x != trace.rend(); ++x) {
                    auto &aln = **x;
                    if (aln.ok) {
                        if (ok_alns == 0) {
                            query_start = aln.j;
                            target_start = aln.i;
                        }
                        ++ok_alns;
                        if (query_end && aln.j > query_end) {
                            const int len = aln.j - query_end;
                            for (uint64_t i = 0; i < len; ++i) {
                                rawv.push_back('I');
                            }
                        }
                        if (target_end && aln.i > target_end) {
                            const int len = aln.i - target_end;
                            for (uint64_t i = 0; i < len; ++i) {
                                rawv.push_back('D');
                            }
                        }
                        uint64_t target_aligned_length = 0;
                        uint64_t query_aligned_length = 0;
                        const int start_idx = aln.edit_cigar.begin_offset;
                        const int end_idx = aln.edit_cigar.end_offset;
                        for (int i = start_idx; i < end_idx; i++) {
                            const auto &c = aln.edit_cigar.cigar_ops[i];
                            switch (c) {
                                case 'M':
                                case 'X':
                                    ++query_aligned_length;
                                    ++target_aligned_length;
                                    break;
                                case 'I':
                                    ++query_aligned_length;
                                    break;
                                case 'D':
                                    ++target_aligned_length;
                                    break;
                                default:
                                    break;
                            }
                            rawv.push_back(c);
                        }
                        query_end = aln.j + query_aligned_length;
                        target_end = aln.i + target_aligned_length;
                    }
                    // clean up
                    delete *x;
                }

#ifdef VALIDATE_WFA_WFLIGN
                if (!validate_trace(rawv, query, target - target_pointer_shift,
                            query_length, target_length,
                            query_start, target_start)) {
                    std::cerr
                        << "cigar failure in rawv (at end) "
                        << "\t" << query_name << "\t" << query_total_length
                        << "\t"
                        << query_offset + (query_is_rev
                                               ? query_length - query_end
                                               : query_start)
                        << "\t"
                        << query_offset + (query_is_rev
                                               ? query_length - query_start
                                               : query_end)
                        << "\t" << (query_is_rev ? "-" : "+") << "\t"
                        << target_name << "\t" << target_total_length << "\t"
                        << target_offset - target_pointer_shift + target_start
                        << "\t" << target_offset + target_end << std::endl;
                    exit(1);
                }
#endif

#ifdef WFLIGN_DEBUG
                std::cerr << "[wflign::wflign_affine_wavefront] eroding "
                             "traceback at k="
                          << erode_k << std::endl;
#endif

                // erode by removing matches < k
                for (uint64_t i = 0; i < rawv.size();) {
                    if (rawv[i] == 'M' || rawv[i] == 'X') {
                        uint64_t j = i;
                        while (++j < rawv.size() &&
                               (rawv[j] == 'M' || rawv[j] == 'X')) {
                        }
                        if (j - i < erode_k) {
                            while (i < j) {
                                erodev.push_back('D');
                                erodev.push_back('I');
                                ++i;
                            }
                        } else {
                            while (i < j) {
                                erodev.push_back(rawv[i++]);
                            }
                        }
                    } else {
                        erodev.push_back(rawv[i++]);
                    }
                }
            }

#ifdef VALIDATE_WFA_WFLIGN
            if (!validate_trace(erodev, query, target - target_pointer_shift,
                        query_length, target_length, query_start,
                        target_start)) {
                std::cerr << "cigar failure in erodev "
                          << "\t" << query_name << "\t" << query_total_length
                          << "\t"
                          << query_offset + (query_is_rev
                                                 ? query_length - query_end
                                                 : query_start)
                          << "\t"
                          << query_offset + (query_is_rev
                                                 ? query_length - query_start
                                                 : query_end)
                          << "\t" << (query_is_rev ? "-" : "+") << "\t"
                          << target_name << "\t" << target_total_length << "\t"
                          << target_offset - target_pointer_shift + target_start
                          << "\t" << target_offset + target_end << std::endl;
                exit(1);
            }
#endif

#ifdef WFLIGN_DEBUG
            std::cerr << "[wflign::wflign_affine_wavefront] normalizing eroded "
                         "traceback"
                      << std::endl;
#endif

            // normalize: sort so that I<D and otherwise leave it as-is
            sort_indels(erodev);

#ifdef WFLIGN_DEBUG
            std::cerr << "[wflign::wflign_affine_wavefront] got normalized "
                 "eroded traceback: ";
            for (auto c : erodev) {
                std::cerr << c;
            }
            std::cerr << std::endl;
#endif

            //std::cerr << "FIRST PATCH ROUND" << std::endl;
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            patching(erodev, tracev, 4096, 8, 512, true);

#ifdef VALIDATE_WFA_WFLIGN
            if (!validate_trace(tracev, query,
                        target - target_pointer_shift, query_length,
                        target_length, query_start, target_start)) {
                std::cerr << "cigar failure in tracev "
                          << "\t" << query_name << "\t" << query_total_length
                          << "\t"
                          << query_offset + (query_is_rev
                                                 ? query_length - query_end
                                                 : query_start)
                          << "\t"
                          << query_offset + (query_is_rev
                                                 ? query_length - query_start
                                                 : query_end)
                          << "\t" << (query_is_rev ? "-" : "+") << "\t"
                          << target_name << "\t" << target_total_length << "\t"
                          << target_offset - target_pointer_shift + target_start
                          << "\t" << target_offset + target_end << std::endl;
                exit(1);
            }
#endif

        }
    }

    // normalize the indels
    sort_indels(tracev);

#ifdef WFLIGN_DEBUG
    std::cerr << "[wflign::wflign_affine_wavefront] got full patched traceback: ";
    for (auto c : tracev) {
        std::cerr << c;
    }
    std::cerr << std::endl;
#endif

#ifdef VALIDATE_WFA_WFLIGN
//            std::cerr << "query_length: " << query_length << std::endl;
//            std::cerr << "target_length: " << target_length <<
//            std::endl; std::cerr << "query_start: " << query_start <<
//            std::endl; std::cerr << "target_start: " << target_start <<
//            std::endl;

    if (!validate_trace(tracev, query, target - target_pointer_shift,
                    query_length, target_length + 2 * wflign_max_len_minor, query_start,
                    target_start)) {
        std::cerr
            << "cigar failure at alignment (before head/tail del trimming) "
            << "\t" << query_name << "\t" << query_total_length << "\t"
            << query_offset +
                   (query_is_rev ? query_length - query_end : query_start)
            << "\t"
            << query_offset +
                   (query_is_rev ? query_length - query_start : query_end)
            << "\t" << (query_is_rev ? "-" : "+") << "\t" << target_name << "\t"
            << target_total_length << "\t"
            << target_offset - target_pointer_shift + target_start << "\t"
            << target_offset + target_end << std::endl;
        exit(1);
    }
#endif

    // trim deletions at start and end of tracev
    uint64_t begin_offset;
    uint64_t end_offset;
    {
        uint64_t trim_del_first;
        uint64_t trim_del_last;

        // 1.) sort initial ins/del to put del < ins
        auto first_non_indel = tracev.begin();
        while (first_non_indel != tracev.end() &&
               (*first_non_indel == 'D' || *first_non_indel == 'I')) {
            ++first_non_indel;
        }
        std::sort(tracev.begin(), first_non_indel,
                  [](char a, char b) { return a < b; });
        // 2.) find first non-D in tracev --> tracev_begin
        //   a.) add to target_start this count
        auto first_non_del = tracev.begin();
        while (first_non_del != tracev.end() && *first_non_del == 'D') {
            ++first_non_del;
        }
        trim_del_first = std::distance(tracev.begin(), first_non_del);
        target_start += trim_del_first;

        // 3.) count D's at end of tracev --> tracev_end
        //   b.) subtract from target_end this count
        auto last_non_del = tracev.rbegin();
        while (last_non_del != tracev.rend() && *last_non_del == 'D') {
            ++last_non_del;
        }
        trim_del_last = std::distance(tracev.rbegin(), last_non_del);
        target_end -= trim_del_last;

        begin_offset = trim_del_first;
        end_offset = tracev.size() - trim_del_last;
    }

    /*
#ifdef VALIDATE_WFA_WFLIGN
    if (!validate_trace(tracev, query, target - target_pointer_shift,
query_length, target_length, query_start, target_start)) { std::cerr <<
"cigar failure at alignment (after head/tail del trimming) "
                  << "\t" << query_name
                  << "\t" << query_total_length
                  << "\t" << query_offset + (query_is_rev ? query_length -
query_end : query_start)
                  << "\t" << query_offset + (query_is_rev ? query_length -
query_start : query_end)
                  << "\t" << (query_is_rev ? "-" : "+")
                  << "\t" << target_name
                  << "\t" << target_total_length
                  << "\t" << target_offset - target_pointer_shift + target_start
                  << "\t" << target_offset + target_end << std::endl;
        exit(1);
    }
#endif
    */

#ifdef WFA_PNG_TSV_TIMING
    bool emit_png = !prefix_wavefront_plot_in_png->empty() && wfplot_max_size > 0;
    if (emit_png) {
        const int pattern_length = (int)query_length;
        const int text_length = (int)target_length;

        const int wfplot_vmin = 0, wfplot_vmax = pattern_length;
        const int wfplot_hmin = 0, wfplot_hmax = text_length;

        int v_max = wfplot_vmax - wfplot_vmin;
        int h_max = wfplot_hmax - wfplot_hmin;

        const algorithms::color_t COLOR_MASH_MISMATCH = { 0xffefefef };
        const algorithms::color_t COLOR_WFA_MISMATCH = { 0xffff0000 };
        const algorithms::color_t COLOR_WFA_MATCH = { 0xff00ff00 };

        const double scale = std::min(1.0, (double)wfplot_max_size / (double)std::max(v_max, h_max));

        const uint64_t width = (int)(scale * (double)v_max);
        const uint64_t height = (int)(scale * (double)h_max);
        const double source_width = (double)width;
        const double source_height = (double)height;

        const double x_off = 0, y_off = 0;
        const double line_width = 1.0;
        const double source_min_x = 0, source_min_y = 0;

        auto plot_point = (v_max <= wfplot_max_size && h_max <= wfplot_max_size)
                          ? [](const algorithms::xy_d_t &point, algorithms::atomic_image_buf_t& image, const algorithms::color_t &color) {
                    image.layer_pixel(point.x, point.y, color);
                }
                          : [](const algorithms::xy_d_t &point, algorithms::atomic_image_buf_t& image, const algorithms::color_t &color) {
                    wflign::algorithms::wu_calc_wide_line(
                            point, point,
                            color,
                            image);
                };

        // Plot the traceback
        {
            algorithms::atomic_image_buf_t image(width, height,
                                                 source_width, source_height,
                                                 source_min_x, source_min_y);


            uint64_t v = query_start; // position in the pattern
            uint64_t h = target_start; // position in the text
            int64_t last_v = -1;
            int64_t last_h = -1;
            for (const auto& c : tracev) {
                switch (c) {
                    case 'M':
                    case 'X':
                        ++v;
                        ++h;
                        {
                            uint64_t _v = v;
                            uint64_t _h = h;
                            if ((_v != last_v && _h != last_h)
                                && _v >= wfplot_vmin && _v <= wfplot_vmax
                                && _h >= wfplot_hmin && _h <= wfplot_hmax) {
                                algorithms::xy_d_t xy0 = {
                                        (_v * scale) - x_off,
                                        (_h * scale) + y_off
                                };
                                xy0.into(source_min_x, source_min_y,
                                         source_width, source_height,
                                         0, 0,
                                         width, height);
                                plot_point(xy0, image, COLOR_WFA_MATCH);
                                last_v = _v;
                                last_h = _h;
                            }
                        }
                        break;
                    case 'I':
                        ++v;
                        break;
                    case 'D':
                        ++h;
                        break;
                    default:
                        break;
                }
                //std::cerr << "plot cell " << v << "," << h << std::endl;
            }

            auto bytes = image.to_bytes();
            const std::string filename = *prefix_wavefront_plot_in_png +
                                         query_name + "_" + std::to_string(query_offset) + "_" + std::to_string(query_offset+query_length) + " _ " + (query_is_rev ? "-" : "+") +
                                         "_" + target_name + "_" + std::to_string(target_offset) + "_" + std::to_string(target_offset+target_length) + ".2.trace.png";
            encodeOneStep(filename.c_str(), bytes, width, height);
        }
    }
#endif

    // convert trace to cigar, get correct start and end coordinates
    char *cigarv = alignment_to_cigar(
            tracev, begin_offset, end_offset,
            total_target_aligned_length, total_query_aligned_length, matches,
            mismatches, insertions, inserted_bp, deletions, deleted_bp);

    // Calculate metrics from the compressed CIGAR string
    uint64_t cigar_matches = 0;
    uint64_t cigar_mismatches = 0;
    uint64_t cigar_insertions = 0;
    uint64_t cigar_inserted_bp = 0;
    uint64_t cigar_deletions = 0;
    uint64_t cigar_deleted_bp = 0;
    uint64_t cigar_refAlignedLength = 0;
    uint64_t cigar_qAlignedLength = 0;

    process_compressed_cigar(
        std::string(cigarv),
        cigar_matches,
        cigar_mismatches,
        cigar_insertions,
        cigar_inserted_bp,
        cigar_deletions,
        cigar_deleted_bp,
        cigar_refAlignedLength,
        cigar_qAlignedLength);

    const double gap_compressed_identity =
            (double)cigar_matches /
            (double)(cigar_matches + cigar_mismatches + cigar_insertions + cigar_deletions);

    const uint64_t edit_distance = cigar_mismatches + cigar_inserted_bp + cigar_deleted_bp;

    // identity over the full block
    const double block_identity =
            (double)cigar_matches / (double)(cigar_matches + edit_distance);

    // Apply filtering using configurable thresholds
    if (gap_compressed_identity >= min_identity && cigar_qAlignedLength >= min_alignment_length && block_identity >= min_block_identity) {
        auto write_tag_and_md_string = [&](std::ostream &out, const char *c,
                                           const int target_start) {
            out << "MD:Z:";

            char last_op = '\0';
            int last_len = 0;
            int t_off = target_start, l_MD = 0;
            int l = 0;
            int x = 0;
            while (c[x] != '\0') {
                while (isdigit(c[x]))
                    ++x;
                char op = c[x];
                int len = 0;
                std::from_chars(c + l, c + x, len);
                l = ++x;
                if (last_len) {
                    if (last_op == op) {
                        len += last_len;
                    } else {
                        // std::cerr << t_off << "   " << last_len << last_op <<
                        // std::endl;

                        if (last_op == '=') {
                            l_MD += last_len;
                            t_off += last_len;
                        } else if (last_op == 'X') {
                            for (uint64_t ii = 0; ii < last_len; ++ii) {
                                out << l_MD
                                    << target[t_off + ii - target_pointer_shift];
                                l_MD = 0;
                            }

                            t_off += last_len;
                        } else if (last_op == 'D') {
                            out << l_MD << "^";
                            for (uint64_t ii = 0; ii < last_len; ++ii) {
                                out << target[t_off + ii - target_pointer_shift];
                            }

                            l_MD = 0;
                            t_off += last_len;
                        }
                    }
                }
                last_op = op;
                last_len = len;
            }

            if (last_len) {
                // std::cerr << t_off << "   " << last_len << last_op << std::endl;

                if (last_op == '=') {
                    out << last_len + l_MD;
                } else if (last_op == 'X') {
                    for (uint64_t ii = 0; ii < last_len; ++ii) {
                        out << l_MD << target[t_off + ii - target_pointer_shift];
                        l_MD = 0;
                    }
                    out << "0";
                } else if (last_op == 'I') {
                    out << l_MD;
                } else if (last_op == 'D') {
                    out << l_MD << "^";
                    for (uint64_t ii = 0; ii < last_len; ++ii) {
                        out << target[t_off + ii - target_pointer_shift];
                    }
                    out << "0";
                }
            }
        };

#ifdef WFA_PNG_TSV_TIMING
        const long elapsed_time_patching_ms =
                std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::steady_clock::now() - start_time)
                        .count();

        const std::string timings_and_num_alignements =
                "wt:i:" + std::to_string(elapsed_time_wflambda_ms) +
                "\tpt:i:" + std::to_string(elapsed_time_patching_ms) +
                "\taa:i:" + std::to_string(num_alignments) +
                "\tap:i:" + std::to_string(num_alignments_performed);
#endif

        if (paf_format_else_sam) {
            out << query_name << "\t" << query_total_length << "\t"
                << query_offset +
                   (query_is_rev ? query_length - query_end : query_start)
                << "\t"
                << query_offset +
                   (query_is_rev ? query_length - query_start : query_end)
                << "\t" << (query_is_rev ? "-" : "+") << "\t" << target_name
                << "\t" << target_total_length << "\t"
                << target_offset - target_pointer_shift + target_start << "\t"
                << target_offset + target_end << "\t" << cigar_matches << "\t"
                << cigar_matches + cigar_mismatches + cigar_inserted_bp + cigar_deleted_bp
                << "\t"
                << std::round(float2phred(1.0 - block_identity))
                //<< "\t" << "as:i:" << total_score
                << "\t"
                << "gi:f:" << gap_compressed_identity << "\t"
                << "bi:f:"
                << block_identity
                //<< "\t" << "md:f:" << mash_dist_sum / trace.size()
                //<< "\t" << "ma:i:" << matches
                //<< "\t" << "mm:i:" << mismatches
                //<< "\t" << "ni:i:" << insertions
                //<< "\t" << "ii:i:" << inserted_bp
                //<< "\t" << "nd:i:" << deletions
                //<< "\t" << "dd:i:" << deleted_bp
                << "\t"
                << "md:f:" << mashmap_estimated_identity;

            if (emit_md_tag) {
                out << "\t";

                write_tag_and_md_string(out, cigarv, target_start);
            }

#ifdef WFA_PNG_TSV_TIMING
            out << "\t" << timings_and_num_alignements << "\t"
                << "cg:Z:" << cigarv << "\n";
#else
            out << "\t" << "cg:Z:" << cigarv << "\n";
#endif
        } else {
            out << "@PG\tID:wfmash\tPN:wfmash\tVN:" << WFLIGN_GIT_VERSION << "\tCL:wfmash\n"
                << query_name                          // Query template NAME
                << "\t" << (query_is_rev ? "16" : "0") // bitwise FLAG
                << "\t" << target_name // Reference sequence NAME
                << "\t"
                << target_offset - target_pointer_shift + target_start +
                   1 // 1-based leftmost mapping POSition
                << "\t"
                << std::round(
                        float2phred(1.0 - block_identity)) // MAPping Quality
                << "\t";

            // CIGAR
            const uint64_t query_start_pos =
                    query_offset +
                    (query_is_rev ? query_length - query_end : query_start);
            const uint64_t query_end_pos =
                    query_offset +
                    (query_is_rev ? query_length - query_start : query_end);

            if (query_start_pos > 0) {
                out << query_start_pos << "H";
            }
            out << cigarv;
            if (query_total_length > query_end_pos) {
                out << (query_total_length - query_end_pos) << "H";
            }

            out << "\t"
                << "*" // Reference name of the mate/next read
                << "\t"
                << "0" // Position of the mate/next read
                << "\t"
                << "0" // observed Template LENgth
                << "\t";

            // segment SEQuence
            if (no_seq_in_sam) {
                out << "*";
            } else {
                for (uint64_t p = query_start; p < query_end; ++p) {
                    out << query[p];
                }
            }

            out << "\t"
                << "*" // ASCII of Phred-scaled base QUALity+33
                << "\t"
                << "NM:i:"
                << edit_distance
                //<< "\t" << "AS:i:" << total_score
                << "\t"
                << "gi:f:" << gap_compressed_identity << "\t"
                << "bi:f:"
                << block_identity
                //<< "\t" << "md:f:" << mash_dist_sum / trace.size()
                //<< "\t" << "ma:i:" << matches
                //<< "\t" << "mm:i:" << mismatches
                //<< "\t" << "ni:i:" << insertions
                //<< "\t" << "ii:i:" << inserted_bp
                //<< "\t" << "nd:i:" << deletions
                //<< "\t" << "dd:i:" << deleted_bp
                << "";

            if (emit_md_tag) {
                out << "\t";

                write_tag_and_md_string(out, cigarv, target_start);
            }
#ifdef WFA_PNG_TSV_TIMING
            out << "\t" << timings_and_num_alignements << "\n";
#else
            out << "\n";
#endif
        }
    }

    // Clean up
    free(cigarv);
    
    // Write SAM format alignments and clean up trace
    if (!paf_format_else_sam) {
        // Clean up the trace alignments since we're done with them
        for (auto* aln : trace) {
            if (aln != nullptr) {
                delete aln;
            }
        }
        
        // Write the patch alignments
        for (auto& patch_aln : multi_patch_alns) {
            write_alignment_sam(
                out, patch_aln, "", query_name, query_total_length,
                query_offset, query_length, query_is_rev,
                target_name, target_total_length, target_offset, target_length,
                min_identity, min_alignment_length, min_block_identity, mashmap_estimated_identity,
                no_seq_in_sam, emit_md_tag, query, target, target_pointer_shift,
                0, 0, 0);
        }
        
        // Clean up patch alignments after writing
        for (auto& patch_aln : multi_patch_alns) {
            free(patch_aln.edit_cigar.cigar_ops);
            patch_aln.edit_cigar.cigar_ops = nullptr;
        }
        multi_patch_alns.clear();
    } else {
        // write how many reverse complement alignments were found
        //std::cerr << "got " << rev_patch_alns.size() << " rev patch alns" << std::endl;
        for (auto& patch_aln : multi_patch_alns) {
            // write_alignment_paf only writes anything if aln.ok. We need to guard the manual tag writing below with the same conditional to avoid writing an invalid PAF.
            bool wrote = write_alignment_paf(
                    out,
                    patch_aln,
                    "",
                    query_name,
                    query_total_length,
                    query_offset,
                    query_length,
                    query_is_rev,
                    target_name,
                    target_total_length,
                    target_offset,
                    target_length,
                    min_identity,
                    min_alignment_length,
                    min_block_identity,
                    mashmap_estimated_identity,
                    0, 0, 0,
                    false,  // Don't add an endline after each alignment
                    true);  // This is a reverse complement alignment
            if (wrote) {
                // write tag indicating that this is a multipatch alignment
                out << "\t" << "pt:Z:true" << "\t"
                    // and if the patch is inverted as well
                    << "iv:Z:" << (patch_aln.is_rev ? "true" : "false") << "\n";
            }
        }
    }
    out << std::flush;
}

void write_tag_and_md_string(
    std::ostream &out,
    const char *cigar_ops,
    const int cigar_start,
    const int cigar_end,
    const int target_start,
    const char *target) {

    out << "MD:Z:";

    char last_op = '\0';
    int last_len = 0;
    int t_off = target_start, l_MD = 0;
    int l = cigar_start;
    int x = cigar_start;

    while (x < cigar_end) {
        // Parse the length digits
        while (x < cigar_end && isdigit(cigar_ops[x]))
            ++x;

        char op = cigar_ops[x];
        int len = 0;

        // Convert the substring [l, x) to an integer
        auto result = std::from_chars(cigar_ops + l, cigar_ops + x, len);
        l = ++x;

        // Process the previous operation if there was one
        if (last_len) {
            if (last_op == op) {
                len += last_len;
            } else {
                // Handle the previous operation based on its type
                if (last_op == '=' || last_op == 'M') {
                    l_MD += last_len;
                    t_off += last_len;
                } else if (last_op == 'X') {
                    for (uint64_t ii = 0; ii < last_len; ++ii) {
                        int64_t idx = t_off + ii ;
                        out << l_MD << target[idx];
                        l_MD = 0;
                    }
                    t_off += last_len;
                } else if (last_op == 'D') {
                    out << l_MD << "^";
                    for (uint64_t ii = 0; ii < last_len; ++ii) {
                        int64_t idx = t_off + ii;
                        out << target[idx];
                    }
                    l_MD = 0;
                    t_off += last_len;
                }
            }
        }
        last_op = op;
        last_len = len;
    }

    // Process the last operation
    if (last_len) {
        if (last_op == '=' || last_op == 'M') {
            out << last_len + l_MD;
        } else if (last_op == 'X') {
            for (uint64_t ii = 0; ii < last_len; ++ii) {
                int64_t idx = t_off + ii;
                out << l_MD << target[idx];
                l_MD = 0;
            }
            out << "0";
        } else if (last_op == 'I') {
            out << l_MD;
        } else if (last_op == 'D') {
            out << l_MD << "^";
            for (uint64_t ii = 0; ii < last_len; ++ii) {
                int64_t idx = t_off + ii;
                out << target[idx];
            }
            out << "0";
        }
    }
}

void write_alignment_sam(
    std::ostream &out,
    const alignment_t& patch_aln,
    const std::string& cigar_str,
    const std::string& query_name,
    const uint64_t& query_total_length,
    const uint64_t& query_offset,
    const uint64_t& query_length,
    const bool& query_is_rev,
    const std::string& target_name,
    const uint64_t& target_total_length,
    const uint64_t& target_offset,
    const uint64_t& target_length,
    const float& min_identity,
    const uint64_t& min_alignment_length,
    const float& min_block_identity,
    const float& mashmap_estimated_identity,
    const bool& no_seq_in_sam,
    const bool& emit_md_tag,
    const char* query,
    const char* target,
    const int64_t& target_pointer_shift,
    const int32_t& chain_id,
    const int32_t& chain_length,
    const int32_t& chain_pos
) {

    if (cigar_str == "") { std::cerr << "[wflign_patch] unsupported codepath" << std::endl; exit(1); }

    uint64_t patch_matches = 0;
    uint64_t patch_mismatches = 0;
    uint64_t patch_insertions = 0;
    uint64_t patch_inserted_bp = 0;
    uint64_t patch_deletions = 0;
    uint64_t patch_deleted_bp = 0;
    uint64_t patch_refAlignedLength = 0;
    uint64_t patch_qAlignedLength = 0;

    process_compressed_cigar(
        cigar_str,
        patch_matches,
        patch_mismatches,
        patch_insertions,
        patch_inserted_bp,
        patch_deletions,
        patch_deleted_bp,
        patch_refAlignedLength,
        patch_qAlignedLength);

    // Trim deletions and get new coordinates
    auto [trimmed_cigar, new_ref_start, new_ref_end, new_query_start, new_query_end] = 
        trim_indels(cigar_str, target_offset + patch_aln.i, target_offset + patch_aln.i + patch_refAlignedLength,
                      query_offset + patch_aln.j, query_offset + patch_aln.j + patch_qAlignedLength);

    // Recompute metrics with trimmed CIGAR
    process_compressed_cigar(
        trimmed_cigar,
        patch_matches,
        patch_mismatches,
        patch_insertions,
        patch_inserted_bp,
        patch_deletions,
        patch_deleted_bp,
        patch_refAlignedLength,
        patch_qAlignedLength);

    char* patch_cigar = strdup(trimmed_cigar.c_str());

    double patch_gap_compressed_identity = (double)patch_matches /
        (double)(patch_matches + patch_mismatches + patch_insertions + patch_deletions);
    double patch_block_identity = (double)patch_matches /
        (double)(patch_matches + patch_mismatches + patch_inserted_bp + patch_deleted_bp);

    /*
    const uint64_t query_start_pos =
        (query_is_rev ? query_offset + query_length : query_offset);
    const uint64_t query_end_pos =
        (query_is_rev ? query_offset : query_offset + query_length);
    */

    // Apply filtering using configurable thresholds
    if (patch_gap_compressed_identity >= min_identity && patch_qAlignedLength >= min_alignment_length && patch_block_identity >= min_block_identity) {

        out << query_name << "\t"
            << (query_is_rev ^ patch_aln.is_rev ? "16" : "0") << "\t"
            << target_name << "\t"
            << new_ref_start + 1 << "\t"
            << std::round(float2phred(1.0 - patch_block_identity)) << "\t"
            << patch_cigar << "\t"
            << "*\t0\t0\t";

        if (no_seq_in_sam) {
            out << "*";
        } else {
            std::stringstream seq;
            // The new patch_aln.j is "new_query_start - query_offset"
            for (uint64_t p = (new_query_start - query_offset); p < (new_query_start - query_offset) + patch_qAlignedLength; ++p) {
                seq << query[p];
            }
            if (patch_aln.is_rev) {
                // reverse complement
                out << reverse_complement(seq.str());
            } else {
                out << seq.str();
            }
        }
        out << "\t*\t"
            << "NM:i:" << (patch_mismatches + patch_inserted_bp + patch_deleted_bp) << "\t"
            << "gi:f:" << patch_gap_compressed_identity << "\t"
            << "bi:f:" << patch_block_identity << "\t"
            << "md:f:" << mashmap_estimated_identity
            << (chain_length > 0 ? ("\t" + std::string("ci:i:") + std::to_string(chain_id)) : "")
            << (chain_length > 0 ? ("\t" + std::string("ch:Z:") + std::to_string(chain_id) + "." + std::to_string(chain_length) + "." + std::to_string(chain_pos)) : "");
            //<< "\t" << "pt:Z:true" <<
            //<< "\t" << "iv:Z:" << (patch_aln.is_rev ? "true" : "false");

        if (emit_md_tag) {
            out << "\t";

            // target_start is 0, as we are working on a subset of the target sequence
            write_tag_and_md_string(out, patch_cigar, 0, strlen(patch_cigar), 
                                    0 + patch_aln.i,
                                    target);
        }

        out << "\n";
    }

    free(patch_cigar);
}

bool write_alignment_paf(
        std::ostream& out,
        const alignment_t& aln,
        const std::string& cigar_str,
        const std::string& query_name,
        const uint64_t& query_total_length,
        const uint64_t& query_offset, // query offset on the forward strand
        const uint64_t& query_length, // used to compute the coordinates for reversed alignments
        const bool& query_is_rev, // if the base homology mapping is in the reverse complement orientation
        const std::string& target_name,
        const uint64_t& target_total_length,
        const uint64_t& target_offset,
        const uint64_t& target_length, // unused
        const float& min_identity,
        const uint64_t& min_alignment_length,
        const float& min_block_identity,
        const float& mashmap_estimated_identity,
        const int32_t& chain_id,
        const int32_t& chain_length,
        const int32_t& chain_pos,
        const bool& with_endline,
        const bool& is_rev_patch) {
    bool ret = false;  // return true if we wrote the alignment
    if (cigar_str == "") { std::cerr << "[wflign_patch] unsupported codepath" << std::endl; exit(1); }

    if (aln.ok) {
        uint64_t matches = 0;
        uint64_t mismatches = 0;
        uint64_t insertions = 0;
        uint64_t inserted_bp = 0;
        uint64_t deletions = 0;
        uint64_t deleted_bp = 0;
        uint64_t refAlignedLength = 0;
        uint64_t qAlignedLength = 0;

        process_compressed_cigar(
            cigar_str,
            matches,
            mismatches,
            insertions,
            inserted_bp,
            deletions,
            deleted_bp,
            refAlignedLength,
            qAlignedLength);

        // Trim deletions and get new coordinates
        auto [trimmed_cigar, new_ref_start, new_ref_end, new_query_start, new_query_end] = 
            trim_indels(cigar_str, target_offset + aln.i, target_offset + aln.i + refAlignedLength,
                         query_offset + aln.j, query_offset + aln.j + qAlignedLength);

        // Recompute metrics with trimmed CIGAR
        process_compressed_cigar(
            trimmed_cigar,
            matches,
            mismatches,
            insertions,
            inserted_bp,
            deletions,
            deleted_bp,
            refAlignedLength,
            qAlignedLength);

        char* cigar = strdup(trimmed_cigar.c_str());
        size_t alignmentRefPos = new_ref_start - target_offset;
        double gap_compressed_identity =
                (double)matches /
                (double)(matches + mismatches + insertions + deletions);
        double block_identity =
                (double)matches /
                (double)(matches + mismatches + inserted_bp + deleted_bp);

        // Apply filtering using configurable thresholds
        if (gap_compressed_identity >= min_identity && qAlignedLength >= min_alignment_length && block_identity >= min_block_identity) {
            uint64_t q_start, q_end;
            // The new aln.j is "new_query_start - query_offset"
            if (query_is_rev) {
                q_start = query_offset + (query_length - (new_query_start - query_offset) - qAlignedLength);
                q_end = query_offset + (query_length - (new_query_start - query_offset));
            } else {
                q_start = new_query_start; // query_offset + (new_query_start - query_offset);
                q_end = new_query_start + qAlignedLength; // query_offset + (new_query_start - query_offset) + qAlignedLength;
            }

            out << query_name << "\t" << query_total_length << "\t" << q_start
                << "\t" << q_end << "\t"
                << (aln.is_rev ^ query_is_rev ? "-" : "+") << "\t" << target_name << "\t"
                << target_total_length << "\t"
                << target_offset + alignmentRefPos << "\t"
                << target_offset + alignmentRefPos + refAlignedLength << "\t"
                << matches << "\t" << std::max(refAlignedLength, qAlignedLength)
                << "\t" << std::round(float2phred(1.0 - block_identity)) << "\t"
                //<< "as:i:" << aln.score << "\t"
                << "gi:f:" << gap_compressed_identity << "\t"
                << "bi:f:" << block_identity << "\t"
                << "md:f:" << mashmap_estimated_identity << "\t"
                // << (chain_length > 0 ? ("\t" + std::string("ci:i:") + std::to_string(chain_id)) : "")
                << (chain_length > 0 ? (std::string("ch:Z:") + std::to_string(chain_id) + "." + std::to_string(chain_length) + "." + std::to_string(chain_pos) + "\t") : "")
                //<< "\t" << "ma:i:" << matches
                //<< "\t" << "mm:i:" << mismatches
                //<< "\t" << "ni:i:" << insertions
                //<< "\t" << "bi:i:" << inserted_bp
                //<< "\t" << "nd:i:" << deletions
                //<< "\t" << "bd:i:" << deleted_bp
                << "cg:Z:" << cigar << "\t";
            if (with_endline) {
                out << std::endl;
            }
            ret = true;
        }
        free(cigar);
    }
    return ret;
}

double float2phred(const double& prob) {
    if (prob == 1)
        return 255; // guards against "-0"
    double p = -10 * (double)log10(prob);
    if (p < 0 || p > 255) // int overflow guard
        return 255;
    else
        return p;
}

void sort_indels(std::vector<char>& v) {
    auto f = v.begin();
    while (f != v.end()) {
        auto j = f;
        while (j != v.end() && (*j == 'D' || *j == 'I')) {
            ++j;
        }
        if (j != f) {
            std::sort(f, j, [](char a, char b) { return b < a; });
            f = j;
        } else {
            ++f;
        }
    }
}

    } /* namespace wavefront */
} /* namespace wflign */
