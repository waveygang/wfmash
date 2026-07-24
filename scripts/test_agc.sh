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
