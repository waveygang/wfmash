#!/usr/bin/env bash
# Verify that wfmash produces identical output from a FASTA and from the same sequences
# stored in an AGC archive.
#
# Requires: agc and samtools on PATH, and a wfmash built with AGC support
# (point WFMASH at it; defaults to build/bin/wfmash).
#
# Usage: scripts/test_agc.sh <fasta[.gz]> [extra wfmash args...]
set -euo pipefail

fasta="${1:?usage: test_agc.sh <fasta[.gz]> [wfmash args...]}"
shift || true
wfmash="${WFMASH:-build/bin/wfmash}"

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

# Plain FASTA (+ .fai) for the FASTA run and as the source for the AGC archive.
zcat -f "$fasta" > "$tmp/seq.fa"
samtools faidx "$tmp/seq.fa"
agc create -o "$tmp/seq.agc" "$tmp/seq.fa" > /dev/null 2>&1

# -t 1 makes wfmash deterministic, so the sorted PAFs are directly comparable.
"$wfmash" "$tmp/seq.fa"  -t 1 "$@" 2>/dev/null | sort > "$tmp/fa.paf"
"$wfmash" "$tmp/seq.agc" -t 1 "$@" 2>/dev/null | sort > "$tmp/agc.paf"

if diff -q "$tmp/fa.paf" "$tmp/agc.paf" > /dev/null; then
    echo "[test_agc] OK: FASTA and AGC produce identical output ($(wc -l < "$tmp/fa.paf") mappings) for $fasta"
else
    echo "[test_agc] FAIL: FASTA and AGC output differ for $fasta"
    diff "$tmp/fa.paf" "$tmp/agc.paf" | head
    exit 1
fi

# AGC keeps the whole FASTA header line as the contig name, so an archive built from headers
# that carry a description must still expose the sequences under their first token, the name
# the FASTA reader and the PAF use.
sed 's/^\(>[^ ]*\).*/\1 some description here/' "$tmp/seq.fa" > "$tmp/desc.fa"
samtools faidx "$tmp/desc.fa"
agc create -o "$tmp/desc.agc" "$tmp/desc.fa" > /dev/null 2>&1

"$wfmash" "$tmp/desc.fa"  -t 1 "$@" 2>/dev/null | sort > "$tmp/desc_fa.paf"
"$wfmash" "$tmp/desc.agc" -t 1 "$@" 2>/dev/null | sort > "$tmp/desc_agc.paf"

if [ -s "$tmp/desc_fa.paf" ] && diff -q "$tmp/desc_fa.paf" "$tmp/desc_agc.paf" > /dev/null; then
    echo "[test_agc] OK: headers with descriptions give identical output ($(wc -l < "$tmp/desc_fa.paf") mappings)"
else
    echo "[test_agc] FAIL: FASTA and AGC output differ when headers carry a description"
    diff "$tmp/desc_fa.paf" "$tmp/desc_agc.paf" | head
    exit 1
fi

# A multi-sample archive built the canonical AGC way (one sample per file, plain contig
# names) repeats contig names across samples. Those records must stay distinct: each is
# named contig@sample and must carry its own sample's bases. If they collapsed onto one
# name/sequence, the two copies would look identical and be dropped as self-mappings.
name="$(head -n 1 "$tmp/seq.fa" | sed 's/^>//' | cut -f1 -d' ')"
samtools faidx "$tmp/seq.fa" "$name" | tail -n +2 | tr -d '\n' > "$tmp/one.seq"
{ echo ">shared_ctg"; cat "$tmp/one.seq"; echo; } > "$tmp/sampleX.fa"
# sampleY: the same sequence with every 50th base flipped, so the two still align
{ echo ">shared_ctg"
  awk '{ for (i = 1; i <= length($0); i++) { c = substr($0, i, 1);
         printf "%s", (i % 50 == 0 ? (c == "A" ? "T" : "A") : c) } print "" }' "$tmp/one.seq"
} > "$tmp/sampleY.fa"
agc create -o "$tmp/multi.agc" "$tmp/sampleX.fa" "$tmp/sampleY.fa" > /dev/null 2>&1

"$wfmash" "$tmp/multi.agc" -t 1 "$@" 2>/dev/null > "$tmp/multi.paf" || true
if grep -q 'shared_ctg@sampleX' "$tmp/multi.paf" && grep -q 'shared_ctg@sampleY' "$tmp/multi.paf"; then
    echo "[test_agc] OK: duplicate contig names across samples stay distinct (from $name)"
else
    echo "[test_agc] FAIL: expected shared_ctg@sampleX and shared_ctg@sampleY in the output"
    head "$tmp/multi.paf"
    exit 1
fi
