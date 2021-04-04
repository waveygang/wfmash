# wfmash
threads=16

prefix=chr8.18_22kbps
/usr/bin/time -v wfmash chr8.fa $prefix.fa -t $threads -Na -s 5000 >$prefix.sam && samtools view $prefix.s5000.sam -bS -@ $threads | samtools sort -@ $threads >$prefix.s5000.bam && samtools index $prefix.s5000.bam && samtools coverage $prefix.s5000.bam
/usr/bin/time -v wfmash chr8.fa $prefix.fa -t $threads -Na -s 500 >$prefix.s500.sam && samtools view $prefix.s500.sam -bS -@ $threads | samtools sort -@ $threads >$prefix.s500.bam && samtools index $prefix.s500.bam && samtools coverage $prefix.s500.bam

prefix=chr8_1Mbps_0001
/usr/bin/time -v wfmash chr8.fa $prefix.fa -t $threads -Na -s 50000 >$prefix.s50000.sam && samtools view $prefix.s50000.sam -bS -@ $threads | samtools sort -@ $threads >$prefix.s50000.bam && samtools index $prefix.s50000.bam && samtools coverage $prefix.s50000.bam

# minimap2
prefix=chr8.18_22kbps
/usr/bin/time -v minimap2 chr8.fa $prefix.fa -t $threads -ca >$prefix.mp2.sam && samtools view $prefix.mp2.sam -bS -@ $threads | samtools sort -@ $threads >$prefix.mp2.bam && samtools index $prefix.mp2.bam

prefix=chr8_1Mbps_0001
/usr/bin/time -v minimap2 chr8.fa $prefix.fa -t $threads -ca >$prefix.mp2.sam && samtools view $prefix.mp2.sam -bS -@ $threads | samtools sort -@ $threads >$prefix.mp2.bam && samtools index $prefix.mp2.bam
