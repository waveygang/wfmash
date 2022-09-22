# Prepare pairs of short/medium sequences
# - wfmash -m
wfmash cerevisiae.chrV.fa -m -n 6 -p 70 -s 1k -l 0 -t 16 > cerevisiae.chrV.approx.paf

# - generate query.fa and target.fa + query_target.seq
bash scripts/WFLIGNvsWFA/paf+fasta2seq+fasta.sh cerevisiae.chrV.approx.paf cerevisiae.chrV.fa
samtools faidx cerevisiae.chrV.approx.paf.subset.fa

# -- adapt mapping information
cat cerevisiae.chrV.approx.paf | \
    awk -v OFS='\t' '{print($1":"$3"-"$4"("$5")",$4-$3,0,$4-$3,"+",$6":"$8"-"$9"(+)",$9-$8,0,$9-$8,$10,$11,$12,$13)}' \
    > cerevisiae.chrV.approx.subset.paf

# Run alignments
# - wflign: align query.fa and target.fa, keeping their order (-t 1)
wfmash -n 6 -p 70 -s 1k cerevisiae.chrV.approx.paf.subset.fa -i cerevisiae.chrV.approx.subset.paf -t 1 > cerevisiae.chrV.wflign.paf

# - wfa: align query_target.seq, using 
~/git/WFA2-lib/bin/align_benchmark -a gap-affine-wfa -g 0,4,6,1 --wfa-memory-mode ultralow -i cerevisiae.chrV.approx.paf.seq --output-full cerevisiae.chrV.approx.paf.seq.out
python3 scripts/WFLIGNvsWFA/output2paf.py cerevisiae.chrV.approx.paf.seq.out > cerevisiae.chrV.wfa.paf


# Compare alignments
python3 scripts/WFLIGNvsWFA/compare_alignments.py cerevisiae.chrV.wflign.paf cerevisiae.chrV.wfa.paf
