#!/usr/bin/env bash
#
# Regression test for wfmash issue #399:
#   https://github.com/waveygang/wfmash/issues/399
#
# Two small (~19 Mbp) Cryptococcus neoformans genomes (H99 and Bt22) are
# near-identical 1:1 syntenic assemblies. A correct default all-vs-all mapping
# should therefore be strongly diagonal. Version 0.14.0 was clean; 0.24.x was
# not. This script downloads the genomes on demand, runs wfmash, and asserts
# quantitative quality thresholds so the regression cannot silently return.
#
# Usage:
#   scripts/regression-issue399.sh [path/to/wfmash] [map|full]
#
# Environment:
#   ISSUE399_DATA   cache directory (default: test/data/issue399)
#   ISSUE399_THREADS  threads (default: 8)
#   ISSUE399_MIN_DIAG minimum diagonal bp fraction (default: 0.90)
#   ISSUE399_MIN_ID   minimum bp-weighted identity (default: 0.98)
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

WFMASH="${1:-$REPO_ROOT/build/bin/wfmash}"
MODE="${2:-map}"
THREADS="${ISSUE399_THREADS:-8}"
DATA_DIR="${ISSUE399_DATA:-$REPO_ROOT/test/data/issue399}"
MIN_DIAG="${ISSUE399_MIN_DIAG:-0.90}"
MIN_ID="${ISSUE399_MIN_ID:-0.98}"

H99_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/245/GCF_000149245.1_CNA3/GCF_000149245.1_CNA3_genomic.fna.gz"
BT22_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/056/621/495/GCA_056621495.1_ASM5662149v1/GCA_056621495.1_ASM5662149v1_genomic.fna.gz"
H99="$DATA_DIR/GCF_000149245.1_CNA3_genomic.fna"
BT22="$DATA_DIR/GCA_056621495.1_ASM5662149v1_genomic.fna"

if [[ ! -x "$WFMASH" ]]; then
    echo "ERROR: wfmash binary '$WFMASH' is not executable." >&2
    echo "Build it first (cmake --build build) or pass a path." >&2
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

MODE_FLAG=""
if [[ "$MODE" == "map" ]]; then
    MODE_FLAG="-m"
elif [[ "$MODE" != "full" ]]; then
    echo "ERROR: mode must be 'map' or 'full' (got '$MODE')" >&2
    exit 2
fi

echo "[issue399] running: $WFMASH $MODE_FLAG -t $THREADS <H99> <Bt22>" >&2
# shellcheck disable=SC2086
"$WFMASH" $MODE_FLAG -t "$THREADS" "$H99" "$BT22" >"$PAF" 2>"$LOG"

if [[ ! -s "$PAF" ]]; then
    echo "[issue399] FAIL: empty PAF output (see $LOG)" >&2
    exit 1
fi

# Compute:
#   n         = number of mapping records
#   bp        = total aligned bases
#   diag_frac = fraction of aligned bases assigned to each query chromosome's
#               dominant reference chromosome (1.0 == perfectly syntenic)
#   w_ident   = bp-weighted nucleotide identity (id:f: or gi:f: tag)
read -r N BP DIAG_FRAC W_IDENT < <(
awk -F'\t' '
function alen() {
    a = $11 + 0
    if (a == 0) { a = ($4 + 0) - ($3 + 0); b = ($9 + 0) - ($8 + 0); if (b > a) a = b }
    return a
}
{
    a = alen()
    if (a <= 0) next
    q = $1; r = $6
    bp[q SUBSEP r] += a
    total += a
    n++
    ident = ""
    for (i = 13; i <= NF; i++) {
        if ($i ~ /^gi:f:/) { ident = substr($i, 6); break }
        if ($i ~ /^id:f:/) { ident = substr($i, 6) }
    }
    if (ident != "") { idw += ident * a; idbp += a }
}
END {
    if (total == 0) { print 0, 0, 0, 0; exit }
    for (k in bp) {
        split(k, p, SUBSEP)
        if (bp[k] > best[p[1]]) best[p[1]] = bp[k]
    }
    for (q in best) diag += best[q]
    wi = (idbp > 0) ? idw / idbp : 0
    printf "%d %d %.4f %.4f\n", n, total, diag / total, wi
}
' "$PAF"
)

echo "[issue399] mode=$MODE n=$N bp=$BP diag_frac=$DIAG_FRAC w_ident=$W_IDENT"

# Chain-integrity invariant (issue #399 follow-up): a chained mapping must
# never belong to more than one (query, reference) chromosome pair. A non-zero
# count here means the parallel chain-info array drifted out of sync with the
# mappings after filtering.
CROSS_REF=$(awk -F'\t' '
{
    cid = ""
    for (i = 1; i <= NF; i++)
        if ($i ~ /^ch:Z:/) { cid = substr($i, 6); sub(/\..*/, "", cid) }
    if (cid != "") seen[$1 SUBSEP cid SUBSEP $6] = 1
}
END {
    for (k in seen) { split(k, p, SUBSEP); refs[p[1] SUBSEP p[2]]++ }
    cross = 0
    for (k in refs) if (refs[k] > 1) cross++
    print cross + 0
}
' "$PAF")
echo "[issue399] cross-reference chains=$CROSS_REF (expected 0)"

FAIL=0
if [[ "$CROSS_REF" != "0" ]]; then
    echo "[issue399] FAIL: $CROSS_REF chain(s) span multiple reference chromosomes" >&2
    FAIL=1
fi
if ! awk -v a="$DIAG_FRAC" -v b="$MIN_DIAG" 'BEGIN { exit !(a >= b) }'; then
    echo "[issue399] FAIL: diagonal bp fraction $DIAG_FRAC < $MIN_DIAG" >&2
    FAIL=1
fi
if [[ "$MODE" == "full" ]] || [[ "$MODE" == "map" ]]; then
    if ! awk -v a="$W_IDENT" -v b="$MIN_ID" 'BEGIN { exit !(a >= b) }'; then
        echo "[issue399] FAIL: weighted identity $W_IDENT < $MIN_ID" >&2
        FAIL=1
    fi
fi

if [[ "$FAIL" == "1" ]]; then
    echo "[issue399] FAILED (paf: $PAF, log: $LOG)" >&2
    exit 1
fi

echo "[issue399] PASS"
