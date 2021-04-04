# wfmash
/usr/bin/time -v wfmash chr8_0001.ref chr8_0001.fastq -t 16 -Na > chr8_0001.sam && samtools view chr8_0001.sam -bS | samtools sort > chr8_0001.bam && samtools index chr8_0001.bam
/usr/bin/time -v wfmash chr8_0001.ref chr8_0001.fastq -t 16 -Na -s 1000 > chr8_0001.s1000.sam && samtools view chr8_0001.s1000.sam -bS | samtools sort > chr8_0001.s1000.bam && samtools index chr8_0001.s1000.bam
/usr/bin/time -v wfmash chr8_0001.ref chr8_0001.fastq -t 16 -Na -s 500 > chr8_0001.s500.sam && samtools view chr8_0001.s500.sam -bS | samtools sort > chr8_0001.s500.bam && samtools index chr8_0001.s500.bam

# minimap2
/usr/bin/time -v minimap2 chr8_0001.ref chr8_0001.fastq -t 16 -ca > chr8_0001.mp2.sam && samtools view chr8_0001.mp2.sam -bS | samtools sort > chr8_0001.mp2.bam && samtools index chr8_0001.mp2.bam
