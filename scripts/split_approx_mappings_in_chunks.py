#!/usr/bin/python

# Usage
# Run:
#   python3 split_approx_mappings_in_chunks.py approximate_mappings.paf 4
# It will generate the following files:
#   approximate_mappings.paf.chunk_0.paf
#   approximate_mappings.paf.chunk_1.paf
#   approximate_mappings.paf.chunk_2.paf
#   approximate_mappings.paf.chunk_3.paf

import sys

# The script that takes the approximate mappings, weighs each mapping by computing its length * (1 - estimated identity),
# then creates N new files where the mapping sets have a similar sum of weights.

def split_chunks(l, n):
    result = [[] for i in range(n)]
    sums = [0] * n
    i = 0
    for e in l:
        result[i].append(e)
        sums[i] += e[1]
        i = sums.index(min(sums))
    return result


if __name__ == '__main__':
    path_approx_mappings = sys.argv[1]
    num_of_chunks = int(sys.argv[2])

    rank_to_mapping_dict = {}
    mapping_list = []

    with open(path_approx_mappings) as f:
        for rank, line in enumerate(f):
            # We could avoid keeping everything in memory by reading the file again later
            rank_to_mapping_dict[rank] = line

            _, _, query_start, query_end, _, _, _, target_start, target_end, _, _, _, estimated_identity = line.strip().split('\t')

            num_mapped_bases = max(int(query_end) - int(query_start), int(target_end) - int(target_start))
            estimated_identity = float(estimated_identity.split('id:f:')[1]) / 100.0

            # High divergence makes alignment more difficult
            weight = num_mapped_bases * (1 - estimated_identity)

            mapping_list.append((rank, weight))

	# Chunk the tuples by looking at their weigths
    chunk_list = split_chunks(mapping_list, num_of_chunks)

	# Collect the ranks from the tuples to generate balanced chunks
    for num_chunk, element_list in enumerate(chunk_list):
        with open(path_approx_mappings + f'.chunk_{num_chunk}.paf', 'w') as fw:
            for rank, _ in element_list:
                fw.write(rank_to_mapping_dict[rank])
