from cigar import Cigar

import sys

# It assumes alignments are in the same order in both files
path1_paf = sys.argv[1]
path2_paf = sys.argv[2]

with open(path1_paf) as f1, open(path2_paf) as f2: 
    for line1, line2 in zip(f1, f2):
        # Take info
        query_name_1, query_len_1, query_start_1, query_end_1, query_strand_1, target_name_1, target_len_1, target_start_1, target_end_1 = line1.split('\t')[:9]
        query_name_2, query_len_2, query_start_2, query_end_2, query_strand_2, target_name_2, target_len_2, target_start_2, target_end_2 = line2.split('\t')[:9]

        # Do not check 'query_name_1 != query_name_2 or target_name_1 != target_name_2' for simplicity
        if query_len_1 != query_len_2 or target_len_1 != target_len_2:
            print('Different query and/or target. Something is wrong in input file orders.')
            exit(-1)

        # Take CIGAR strings
        cigar_string_1 = line1.strip().split('cg:Z:')[-1]
        cigar_string_2 = line2.strip().split('cg:Z:')[-1]    

        if not cigar_string_1 or not cigar_string_2:
            print(' Missing CIGAR string(s).')
            exit(-1)

        query_len_1, query_start_1, query_end_1, target_len_1, target_start_1, target_end_1 = [int(x) for x in (query_len_1, query_start_1, query_end_1, target_len_1, target_start_1, target_end_1)]
        query_len_2, query_start_2, query_end_2, target_len_2, target_start_2, target_end_2 = [int(x) for x in (query_len_2, query_start_2, query_end_2, target_len_2, target_start_2, target_end_2)]

        cigar_1_list = list(Cigar(cigar_string_1).items())
        cigar_2_list = list(Cigar(cigar_string_2).items())

        # CIGAR string padding
        if query_start_1 > 0:
            cigar_1_list = [(query_start_1, 'I')] + cigar_1_list
        if query_start_2 > 0:
            cigar_2_list = [(query_start_2, 'I')] + cigar_2_list
        if query_len_1 - query_end_1 > 0:
            cigar_1_list = cigar_1_list + [(query_len_1 - query_end_1, 'I')]
        if query_len_2 - query_end_2 > 0:
            cigar_2_list = cigar_2_list + [(query_len_2 - query_end_2, 'I')]
        if target_start_1 > 0:
            cigar_1_list = [(target_start_1, 'D')] + cigar_1_list
        if target_start_2 > 0:
            cigar_2_list = [(target_start_2, 'D')] + cigar_2_list
        if target_len_1 - target_end_1 > 0:
            cigar_1_list = cigar_1_list + [(target_len_1 - target_end_1, 'D')]
        if target_len_2 - target_end_2 > 0:
            cigar_2_list = cigar_2_list + [(target_len_2 - target_end_2, 'D')]

        # print(cigar_1_list[:10], ' --- ', cigar_1_list[-10:])
        # print(cigar_2_list[:10], ' --- ', cigar_2_list[-10:])

        alignment_pair_set_list = []

        for cigar_x_list in [cigar_1_list, cigar_2_list]:
            alignment_pair_set = set()
            index_q = -1
            index_t = -1
            for (len_op, op) in cigar_x_list:
                if op == 'I':
                    for q in range(index_q, index_q + len_op):
                        alignment_pair_set.add((q + 1, index_t))
                    index_q += len_op
                elif op == 'D':
                    for t in range(index_t, index_t + len_op):
                        alignment_pair_set.add((index_q, t + 1))
                    index_t += len_op
                elif op == 'X' or op == '=' or op == 'M':
                    for i in range(0, len_op):
                        alignment_pair_set.add((index_q + i + 1, index_t + i + 1))
                    index_q += len_op
                    index_t += len_op
                else:
                    print(f'Unexpected CIGAR string operator: {op}.')
                    exit(-1)

            alignment_pair_set_list.append(alignment_pair_set)

        intersection_len = len(alignment_pair_set_list[0].intersection(alignment_pair_set_list[1]))
        union_len = len(alignment_pair_set_list[0].union(alignment_pair_set_list[1]))
        print(float(intersection_len)/float(union_len))
