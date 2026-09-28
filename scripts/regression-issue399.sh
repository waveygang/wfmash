#!/usr/bin/env bash
#
# Regression test for wfmash issue #399:
#   https://github.com/waveygang/wfmash/issues/399
#
# Two small (~19 Mbp) Cryptococcus neoformans genomes (H99 and Bt22) are
# near-identical 1:1 syntenic assemblies. A correct default all-vs-all mapping
# should therefore be strongly diagonal. Version 0.14.0 was clean; 0.24.x was
# not. This script downloads the genomes on demand, runs wfmash, and asserts
# both precision and recall so the regression cannot silently return.
#
# Usage:
#   scripts/regression-issue399.sh [path/to/wfmash] [map|full]
#
# Environment:
#   ISSUE399_DATA     cache directory (default: test/data/issue399)
#   ISSUE399_THREADS  threads (default: 8)
#   ISSUE399_MIN_PRECISION  diag_bp / total_bp (default: 0.90)
#   ISSUE399_MIN_RECALL     retained diag_bp / raw diag_bp (default: 0.95)
#   ISSUE399_MIN_ID         bp-weighted identity (default: 0.98)
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

WFMASH="${1:-$REPO_ROOT/build/bin/wfmash}"
MODE="${2:-map}"
THREADS="${ISSUE399_THREADS:-8}"
DATA_DIR="${ISSUE399_DATA:-$REPO_ROOT/test/data/issue399}"
MIN_PRECISION="${ISSUE399_MIN_PRECISION:-0.90}"
MIN_RECALL="${ISSUE399_MIN_RECALL:-0.95}"
MIN_ID="${ISSUE399_MIN_ID:-0.98}"
STATS="$REPO_ROOT/scripts/paf_stats.py"

H99_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/245/GCF_000149245.1_CNA3/GCF_000149245.1_CNA3_genomic.fna.gz"
BT22_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/056/621/495/GCA_056621495.1_ASM5662149v1/GCA_056621495.1_ASM5662149v1_genomic.fna.gz"
H99="$DATA_DIR/GCF_000149245.1_CNA3_genomic.fna"
BT22="$DATA_DIR/GCA_056621495.1_ASM5662149v1_genomic.fna"

if [[ ! -x "$WFMASH" ]]; then
    echo "ERROR: wfmash binary '$WFMASH' is not executable." >&2
    echo "Build it first (cmake --build build) or pass a path." >&2
    exit 2
fi
if ! command -v python3 >/dev/null 2>&1; then
    echo "ERROR: python3 is required by $STATS" >&2
    exit 2
fi
if [[ "$MODE" != "map" && "$MODE" != "full" ]]; then
    echo "ERROR: mode must be 'map' or 'full' (got '$MODE')" >&2
    exit 2
fi

mkdir -p "$DATA_DIR"

fetch_fasta() {
    local url="$1" out="$2"
    if [[ -s "$out" ]]; then
        return 0
    fi
    echo "[issue399] downloading $(basename "$out") ..." >&2
    if ! curl -fsSL --retry 3 --retry-delay 2 -o "$out.gz" "$url"; then
        echo "[issue399] SKIP: could not download $url (offline?)" >&2
        rm -f "$out.gz"
        return 3
    fi
    gunzip -f "$out.gz"
    if [[ ! -s "$out" ]]; then
        echo "[issue399] ERROR: downloaded FASTA is empty: $out" >&2
        exit 4
    fi
}

if ! fetch_fasta "$H99_URL" "$H99"; then exit 0; fi
if ! fetch_fasta "$BT22_URL" "$BT22"; then exit 0; fi

# wfmash requires FASTA index files (.fai).
index_if_needed() {
    local fa="$1"
    if [[ -s "$fa.fai" ]]; then
        return 0
    fi
    if command -v samtools >/dev/null 2>&1; then
        samtools faidx "$fa"
    else
        echo "[issue399] ERROR: '$fa.fai' missing and samtools not found to create it." >&2
        return 5
    fi
}
index_if_needed "$H99"
index_if_needed "$BT22"

OUT_DIR="$DATA_DIR/out"
mkdir -p "$OUT_DIR"
PAF="$OUT_DIR/issue399.${MODE}.paf"
LOG="$OUT_DIR/issue399.${MODE}.log"
RAW_PAF="$OUT_DIR/issue399.raw.map.paf"
RAW_LOG="$OUT_DIR/issue399.raw.map.log"

MODE_FLAG=""
[[ "$MODE" == "map" ]] && MODE_FLAG="-m"

# Unfiltered mapping run defines the recall ceiling for this genome pair
# (-j 0 disables scaffold filtering; -l 0 disables the block-length filter).
if [[ ! -s "$RAW_PAF" ]]; then
    echo "[issue399] establishing raw ceiling (-m -l 0 -j 0)" >&2
    "$WFMASH" -m -l 0 -j 0 -t "$THREADS" "$H99" "$BT22" >"$RAW_PAF" 2>"$RAW_LOG"
fi

echo "[issue399] running: $WFMASH $MODE_FLAG -t $THREADS <H99> <Bt22>" >&2
# shellcheck disable=SC2086
"$WFMASH" $MODE_FLAG -t "$THREADS" "$H99" "$BT22" >"$PAF" 2>"$LOG"

if [[ ! -s "$PAF" ]] || [[ ! -s "$RAW_PAF" ]]; then
    echo "[issue399] FAIL: empty PAF output (see $LOG)" >&2
    exit 1
fi

stats_value() { sed -n "s/^$2=//p" <<<"$1"; }

RAW_STATS="$(python3 "$STATS" "$RAW_PAF")"
TEST_STATS="$(python3 "$STATS" "$PAF")"

raw_diag_bp="$(stats_value "$RAW_STATS" diag_bp)"
n="$(stats_value "$TEST_STATS" n)"
total_bp="$(stats_value "$TEST_STATS" total_bp)"
diag_bp="$(stats_value "$TEST_STATS" diag_bp)"
precision="$(stats_value "$TEST_STATS" precision)"
diag_cov_frac="$(stats_value "$TEST_STATS" diag_cov_frac)"
w_ident="$(stats_value "$TEST_STATS" w_ident)"
cross_ref="$(stats_value "$TEST_STATS" cross_ref_chains)"

recall="$(awk -v a="$diag_bp" -v b="$raw_diag_bp" 'BEGIN { printf "%.4f", (b > 0 ? a / b : 0) }')"

echo "[issue399] mode=$MODE n=$n total_bp=$total_bp diag_bp=$diag_bp"
echo "[issue399] precision=$precision (min $MIN_PRECISION)  recall_vs_raw=$recall (min $MIN_RECALL)"
echo "[issue399] diag_cov_frac=$diag_cov_frac  w_ident=$w_ident  cross_ref_chains=$cross_ref"

FAIL=0
check_ge() {
    local name="$1" got="$2" want="$3"
    if ! awk -v a="$got" -v b="$want" 'BEGIN { exit !(a >= b) }'; then
        echo "[issue399] FAIL: $name $got < $want" >&2
        FAIL=1
    fi
}
check_ge precision "$precision" "$MIN_PRECISION"
check_ge recall_vs_raw "$recall" "$MIN_RECALL"
check_ge w_ident "$w_ident" "$MIN_ID"

if [[ "$cross_ref" != "0" ]]; then
    echo "[issue399] FAIL: $cross_ref chain(s) span multiple reference chromosomes" >&2
    FAIL=1
fi

if [[ "$FAIL" == "1" ]]; then
    echo "[issue399] FAILED (paf: $PAF, log: $LOG)" >&2
    exit 1
fi

echo "[issue399] PASS"
