#!/bin/bash

# Function to display usage
usage() {
    echo "Usage: $0 -p <paf_file> -f <fasta_file> -o <output_prefix>"
    echo "  -p, --paf       PAF file"
    echo "  -f, --fasta     FASTA file"
    echo "  -o, --output    Output prefix for the new FASTA and adjusted PAF files"
    exit 1
}

# Parse command-line arguments
PARSED_ARGUMENTS=$(getopt -a -n "$0" -o p:f:o: --long paf:,fasta:,output: -- "$@")
VALID_ARGUMENTS=$?
if [ "$VALID_ARGUMENTS" != "0" ]; then
    usage
fi

eval set -- "$PARSED_ARGUMENTS"
while :
do
    case "$1" in
        -p | --paf) PAF_FILE="$2" ; shift 2 ;;
        -f | --fasta) FASTA_FILE="$2" ; shift 2 ;;
        -o | --output) OUTPUT_PREFIX="$2" ; shift 2 ;;
        --) shift ; break ;;
        *) usage ;;
    esac
done

# Validate arguments
if [ -z "$PAF_FILE" ] || [ -z "$FASTA_FILE" ] || [ -z "$OUTPUT_PREFIX" ]; then
    usage
fi

OUTPUT_FASTA="${OUTPUT_PREFIX}.fa"
OUTPUT_PAF="${OUTPUT_PREFIX}.paf"

# Extract the first alignment from the PAF file
read -r line < "$PAF_FILE"

# Parse the alignment fields
IFS=$'\t' read -r -a fields <<< "$line"

# Extract query and target information
query_name="${fields[0]}"
query_start="${fields[2]}"
query_end="${fields[3]}"
target_name="${fields[5]}"
target_start="${fields[7]}"
target_end="${fields[8]}"

# Extract sequences using samtools faidx
samtools faidx "$FASTA_FILE" "$query_name:$query_start-$query_end" > "$OUTPUT_FASTA"
samtools faidx "$FASTA_FILE" "$target_name:$target_start-$target_end" >> "$OUTPUT_FASTA"

qname=$(grep ^">" "$OUTPUT_FASTA" | head -n 1 | sed 's/^>//')
tname=$(grep ^">" "$OUTPUT_FASTA" | tail -n 1 | sed 's/^>//')

# Adjust the PAF file
awk -v qname=$qname \
    -v tname=$tname \
    '{
        if (NF > 0) {
            query_length = $4 - $3;
            target_length = $9 - $8;
            $1 = qname;
            $6 = tname;
            $2 = query_length;
            $3 = 0;
            $4 = query_length;
            $7 = target_length;
            $8 = 0;
            $9 = target_length;
            print $0;
        }
    }' "$PAF_FILE" | tr ' ' '\t' > "$OUTPUT_PAF"

