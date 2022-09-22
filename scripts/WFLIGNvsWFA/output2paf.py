import sys

path_full_output = sys.argv[1]

with open(path_full_output) as f:
    for line in f:
        len_target, len_query, score, target, query, cigar = line.strip().split('\t')

        print('\t'.join(['query', len_query, '0', len_query, '+', 'target', len_target, '0', len_target, 'cg:Z:' + cigar]))
