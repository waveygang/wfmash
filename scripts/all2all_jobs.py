import argparse
import gzip
import itertools

def parse_fasta_index(fasta_file):
    fai_file = fasta_file + '.fai'
    sequences = []
    with open(fai_file, 'r') as file:
        for line in file:
            sequence_name = line.strip().split('\t')[0]
            sequences.append(sequence_name)
    return sequences

def group_sequences(sequences, grouping):
    grouped_sequences = {}
    for sequence in sequences:
        if '#' in sequence:
            fields = sequence.split('#')
            if grouping in ['g', 'genome']:
                group_key = fields[0]
            elif grouping in ['h', 'haplotype']:
                group_key = '#'.join(fields[:2])
            elif grouping in ['c', 'contig']:
                group_key = sequence
            else:
                raise ValueError(f"Invalid grouping: {grouping}")
        else:
            group_key = sequence

        if group_key not in grouped_sequences:
            grouped_sequences[group_key] = []
        grouped_sequences[group_key].append(sequence)

    return grouped_sequences

def generate_pairings(target_grouped_sequences, query_grouped_sequences, num_queries):
    pairings = []
    target_groups = list(target_grouped_sequences.keys())
    query_groups = list(query_grouped_sequences.keys())
    
    for target_group in target_groups:
        query_pool = [group for group in query_groups if group != target_group]
        
        for query_chunk in itertools.zip_longest(*[iter(query_pool)] * num_queries):
            query_chunk = [q for q in query_chunk if q is not None]
            pairings.append((target_group, query_chunk))
    
    return pairings

def main():
    parser = argparse.ArgumentParser(description='Generate pairings or wfmash command lines for all-to-all alignment using PanSN format.')
    parser.add_argument('fasta_file', help='Path to the FASTA file (can be gzipped)')
    parser.add_argument('-n', '--num-queries', type=int, default=4, help='Number of query groups per target group (default: 4)')
    parser.add_argument('-t', '--target-grouping', choices=['g', 'genome', 'h', 'haplotype', 'c', 'contig'], default='haplotype', help='Grouping level for targets: g/genome, h/haplotype, or c/contig (default: haplotype)')
    parser.add_argument('-q', '--query-grouping', choices=['g', 'genome', 'h', 'haplotype', 'c', 'contig'], default='haplotype', help='Grouping level for queries: g/genome, h/haplotype, or c/contig (default: haplotype)')
    parser.add_argument('-o', '--output', help='Output file to save the pairings or command lines')

    args, wfmash_args = parser.parse_known_args()

    # Parse the FASTA index file
    sequences = parse_fasta_index(args.fasta_file)

    # Group sequences based on the specified grouping levels
    target_grouped_sequences = group_sequences(sequences, args.target_grouping)
    query_grouped_sequences = group_sequences(sequences, args.query_grouping)

    # Generate pairings
    pairings = generate_pairings(target_grouped_sequences, query_grouped_sequences, args.num_queries)

    # Save or print the pairings or command lines
    if wfmash_args:
        wfmash_options = ' '.join(wfmash_args)
        if args.output:
            with open(args.output, 'w') as file:
                for target_group, query_groups in pairings:
                    file.write(f"wfmash {wfmash_options} -T {target_group} -Q {','.join(query_groups)}\n")
        else:
            for target_group, query_groups in pairings:
                print(f"wfmash {wfmash_options} -T {target_group} -Q {','.join(query_groups)}")
    else:
        if args.output:
            with open(args.output, 'w') as file:
                for target_group, query_groups in pairings:
                    file.write(f"{target_group}\t{','.join(query_groups)}\n")
        else:
            for target_group, query_groups in pairings:
                print(f"{target_group}\t{','.join(query_groups)}")

if __name__ == '__main__':
    main()