f=cerevisiae.chrV.fa
f=chr8.hprcy1.longs.fa
f=multisaccha.fa

# Prepare pairs of short/medium sequences
# - wfmash -m
wfmash $f -m -n 6 -p 70 -s 1k -l 0 -t 16 | head -n 1 > $f.approx.paf

# - generate query.fa and target.fa + query_target.seq
bash ~/git/wfmash/scripts/WFLIGNvsWFA/paf+fasta2seq+fasta.sh $f.approx.paf $f
samtools faidx $f.approx.paf.subset.fa

# -- adapt mapping information
cat $f.approx.paf | \
    awk -v OFS='\t' '{print($1":"$3"-"$4"("$5")",$4-$3,0,$4-$3,"+",$6":"$8"-"$9"(+)",$9-$8,0,$9-$8,$10,$11,$12,$13)}' \
    > $f.approx.subset.paf

# Run alignments
# - wflign: align query.fa and target.fa, keeping their order (-t 1)
wfmash -n 6 -p 70 -s 1k $f.approx.paf.subset.fa -i $f.approx.subset.paf -t 1 > $f.wflign.paf

# - wfa: align query_target.seq, using 
~/git/WFA2-lib/bin/align_benchmark -a gap-affine-wfa -g 0,4,6,1 --wfa-memory-mode ultralow -i $f.approx.paf.seq --output-full $f.approx.paf.seq.out
python3 ~/git/wfmash/scripts/WFLIGNvsWFA/output2paf.py $f.approx.paf.seq.out > $f.wfa.paf


# Compare alignments
python3 ~/git/wfmash/scripts/WFLIGNvsWFA/compare_alignments.py $f.wflign.paf $f.wfa.paf

# Plots
ls *.wf*.paf | while read f; do ~/git/pafplot/target/release/pafplot $f -d; done