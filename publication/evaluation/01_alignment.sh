# Exit when any command fails
set -eo pipefail

threads=16

prefix=chr8.18_22kbps
/usr/bin/time -v wfmash chr8.fasta $prefix.fasta -t $threads -Na -s 1000 >$prefix.s1000.sam
bash sam2bam.sh $prefix.s1000.sam -@ $threads

/usr/bin/time -v minimap2 chr8.fasta $prefix.fasta -t $threads -ca -x asm20 >$prefix.mp2.sam
bash sam2bam.sh $prefix.mp2.sam -@ $threads

prefix=chr8_1Mbps
/usr/bin/time -v wfmash chr8.fasta $prefix.fasta -t $threads -Na -s 50000 >$prefix.s50000.sam
bash sam2bam.sh $prefix.s50000.sam -@ $threads

/usr/bin/time -v minimap2 chr8.fasta $prefix.fasta -t $threads -ca -x asm20 >$prefix.mp2.sam
bash sam2bam.sh $prefix.mp2.sam -@ $threads
