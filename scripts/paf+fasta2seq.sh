#!/bin/bash

f=$1

wfmash $f -s 1k -p 98 -n 7 -m > $f.paf

cat $f.paf | awk -v OFS='\t' '{print $1, $3, $4, "", "", $5}' > $f.query.bed
cat $f.paf | awk -v OFS='\t' '{print $6, $8, $9, "", "", "+"}' > $f.target.bed

bedtools getfasta -fi $f -bed $f.query.bed -s | grep '^>' -v | sed 's/^/>/' > $f.query.seqs
bedtools getfasta -fi $f -bed $f.target.bed -s | grep '^>' -v | sed 's/^/</'> $f.target.seqs

paste -d '\n' $f.query.seqs $f.target.seqs > $f.seq
