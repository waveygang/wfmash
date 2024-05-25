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
            if grouping == 'genome':
                group_key = fields[0]
            elif grouping == 'haplotype':
                group_key = '#'.join(fields[:2])
            else:
                raise ValueError(f"Invalid grouping: {grouping}")
        else:
            group_key = sequence

        if group_key not in grouped_sequences:
            grouped_sequences[group_key] = []
        grouped_sequences[group_key].append(sequence)

    return grouped_sequences

def generate_pairings(grouped_sequences, num_queries):
    pairings = []
    groups = list(grouped_sequences.keys())
    num_groups = len(groups)
    
    for i in range(num_groups):
        target_group = groups[i]
        query_groups = [group for j, group in enumerate(groups) if j != i]
        
        for query_chunk in itertools.zip_longest(*[iter(query_groups)] * num_queries):
            query_chunk = [q for q in query_chunk if q is not None]
            pairings.append((target_group, query_chunk))
    
    return pairings

def main():
    parser = argparse.ArgumentParser(description='Generate pairings for all-to-all alignment using PanSN format.')
    parser.add_argument('fasta_file', help='Path to the FASTA file (can be gzipped)')
    parser.add_argument('--num-queries', type=int, default=4, help='Number of query groups per target group (default: 4)')
    parser.add_argument('--grouping', choices=['genome', 'haplotype'], default='haplotype', help='Grouping level: genome or haplotype (default: haplotype)')
    parser.add_argument('--output', help='Output file to save the pairings')

    args = parser.parse_args()

    # Parse the FASTA index file
    sequences = parse_fasta_index(args.fasta_file)

    # Group sequences based on the specified grouping level
    grouped_sequences = group_sequences(sequences, args.grouping)

    # Generate pairings
    pairings = generate_pairings(grouped_sequences, args.num_queries)

    # Save or print the pairings
    if args.output:
        with open(args.output, 'w') as file:
            for target_group, query_groups in pairings:
                file.write(f"{target_group},{','.join(query_groups)}\n")
    else:
        for target_group, query_groups in pairings:
            print(f"{target_group},{','.join(query_groups)}")

if __name__ == '__main__':
    main()