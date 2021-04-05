#!/bin/bash

threads=16

prefix=chr8.18_22kbps
/usr/bin/time -v wfmash chr8.fasta $prefix.fasta -t $threads -Na -s 1000 >$prefix.s1000.sam && samtools view $prefix.s1000.sam -bS -@ $threads | samtools sort -@ $threads >$prefix.s1000.bam && samtools index $prefix.s1000.bam && samtools coverage $prefix.s1000.bam

/usr/bin/time -v minimap2 chr8.fasta $prefix.fasta-t $threads -ca >$prefix.mp2.sam && samtools view $prefix.mp2.sam -bS -@ $threads | samtools sort -@ $threads >$prefix.mp2.bam && samtools index $prefix.mp2.bam

prefix=chr8_1Mbps
/usr/bin/time -v wfmash chr8.fasta $prefix.fasta -t $threads -Na -s 50000 >$prefix.s50000.sam && samtools view $prefix.s50000.sam -bS -@ $threads | samtools sort -@ $threads >$prefix.s50000.bam && samtools index $prefix.s50000.bam && samtools coverage $prefix.s50000.bam

/usr/bin/time -v minimap2 chr8.fasta $prefix.fasta -t $threads -ca >$prefix.mp2.sam && samtools view $prefix.mp2.sam -bS -@ $threads | samtools sort -@ $threads >$prefix.mp2.bam && samtools index $prefix.mp2.bam
