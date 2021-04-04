# Downlaod the T2T assembly
wget -c https://s3.amazonaws.com/nanopore-human-wgs/chm13/assemblies/chm13.draft_v1.0.fasta.gz

# Extract the chr8
gunzip chm13.draft_v1.0.fasta.gz
samtools faidx chm13.draft_v1.0.fasta $(grep chr8 chm13.draft_v1.0.fasta.fai | cut -f 1) > chr8.fa

# Simulate [18000, 22000] kbps reads
~/git/PBSIM-PacBio-Simulator/src/pbsim --depth 20 --data-type CLR --accuracy-mean 0.999 --accuracy-min 0.99 --length-min 18000 --length-mean 20000 --length-max 22000 --model_qc ~/git/PBSIM-PacBio-Simulator/data/model_qc_clr --prefix chr8 chr8.fa
