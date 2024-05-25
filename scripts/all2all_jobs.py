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

    return list(grouped_sequences.values())

def generate_pairings(grouped_sequences, num_queries):
    pairings = []
    for target_group in grouped_sequences:
        for query_chunk in itertools.combinations(grouped_sequences, num_queries):
            query_chunk = [seq for group in query_chunk for seq in group if seq not in target_group]
            pairings.append((target_group, query_chunk))
    return pairings

def main():
    parser = argparse.ArgumentParser(description='Generate pairings for all-to-all alignment using PanSN format.')
    parser.add_argument('fasta_file', help='Path to the FASTA file (can be gzipped)')
    parser.add_argument('--num-queries', type=int, default=5, help='Number of query groups per target group (default: 5)')
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
                file.write(f"Target Group: {', '.join(target_group)}\n")
                file.write(f"Query Groups: {', '.join([seq for group in query_groups for seq in group])}\n")
                file.write("\n")
    else:
        for target_group, query_groups in pairings:
            print(f"Target Group: {', '.join(target_group)}")
            print(f"Query Groups: {', '.join([seq for group in query_groups for seq in group])}")
            print()

if __name__ == '__main__':
    main()
