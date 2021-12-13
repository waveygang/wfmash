#!/bin/bash

# Example:
# bash scripts/fasta+paf2vcf.sh data/LPA.subset.fa.gz 100k 300k 98 14 chm13 '_' 16

# Inputs
FASTA=$1
SEGMENT_LEN=$2
BLOCK_LEN=$3
IDENTITY=$4
NUM_HAPLOTYPES=$5
REFERENCE_PREFIX=$6
SEPARATOR=$7
THREADS=$8

# Paths
PAF=$FASTA.s${SEGMENT_LEN}.l${BLOCK_LEN}.p$IDENTITY.n$HAPLOTYPES.k16.paf
GFA=$PAF.k0.B10M.gfa
VCF=$GFA.vcf.gz

echo "All-vs-all alignment"
wfmash $FASTA $FASTA -X -s ${SEGMENT_LEN} -l ${BLOCK_LEN} -p $IDENTITY -n ${NUM_HAPLOTYPES} -k 16 -t $THREADS > $PAF

echo "Graph induction"
# -k 0` is for getting a lossless representation of the pairwise alignment
seqwish -s $FASTA -p $PAF -g $GFA -k 0 -B 10M -t $THREADS -P

echo "Identify variants"
vg deconstruct -e -a -P $REFERENCE_PREFIX -H $SEPARATOR $GFA -t $THREADS | bgzip -c > $VCF && tabix $VCF
