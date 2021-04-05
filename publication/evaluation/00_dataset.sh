#!/bin/bash

# Download the T2T assembly
wget -c https://s3.amazonaws.com/nanopore-human-wgs/chm13/assemblies/chm13.draft_v1.0.fasta.gz

# Extract the chr8
gunzip chm13.draft_v1.0.fasta.gz
samtools faidx chm13.draft_v1.0.fasta
samtools faidx chm13.draft_v1.0.fasta $(grep chr8 chm13.draft_v1.0.fasta.fai | cut -f 1) > chr8.fa

# HiFi read simulation
# [18000, 22000] kbps reads
~/git/PBSIM-PacBio-Simulator/src/pbsim --depth 20 --data-type CLR --accuracy-mean 0.999 --accuracy-min 0.99 --length-min 18000 --length-mean 20000 --length-max 22000 --model_qc ~/git/PBSIM-PacBio-Simulator/data/model_qc_clr --prefix chr8 chr8.fa
sed -n '1~4s/^@/>/p;2~4p' chr8_0001.fastq > chr8.18_22kbps.fasta
rm chr8_0001.ref chr8_0001.fastq

# 1 Mbps reads
~/git/PBSIM-PacBio-Simulator/src/pbsim --depth 20 --data-type CLR --accuracy-mean 0.999 --accuracy-min 0.99 --length-min 1000000 --length-mean 1000000 --length-max 1000000 --model_qc ~/git/PBSIM-PacBio-Simulator/data/model_qc_clr --prefix chr8_1Mbps chr8.fa
sed -n '1~4s/^@/>/p;2~4p' chr8_1Mbps_0001.fastq > chr8.1Mbps.fasta
rm chr8_1Mbps8_0001.ref chr8_1Mbps_0001.fastq
