#!/bin/bash

# Exit when any command fails
set -eo pipefail

path_sam=$1
threads=$2

prefix=$(basename "$path_sam" .sam)

samtools view "$path_sam" -bS | samtools sort -@ "$threads" >"$prefix".bam && samtools index "$prefix".bam
samtools coverage "$prefix".bam
