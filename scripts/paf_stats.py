#!/usr/bin/env python3
"""Summarize PAF quality for near-identical 1:1 genome comparisons.

Emits key=value lines (used by scripts/regression-issue399.sh):

  n                 number of records
  genome            sum of query lengths
  total_bp          total aligned bp
  diag_bp           aligned bp on each query's dominant reference chromosome
  precision         diag_bp / total_bp
  diag_cov          union of query bases mapped to the dominant target
  diag_cov_frac     diag_cov / genome   (recall of the 1:1 alignment)
  all_cov, all_cov_frac
  w_ident           bp-weighted identity (id:f: / gi:f: / 1-dv:f:)
  cross_ref_chains  chains (ch:Z:) spanning >1 reference chromosome (must be 0)

Usage: paf_stats.py alignment.paf
"""
import sys
import collections


def load(path):
    recs = []
    qlen = {}
    with open(path) as f:
        for line in f:
            t = line.rstrip('\n').split('\t')
            if len(t) < 12:
                continue
            q, r = t[0], t[5]
            qs, qe = int(t[2]), int(t[3])
            alen = int(t[10])
            if alen <= 0:
                alen = max(qe - qs, int(t[8]) - int(t[7]))
            ident = None
            chain = None
            for tag in t[12:]:
                if tag.startswith('id:f:'):
                    ident = float(tag[5:])
                elif tag.startswith('gi:f:') and ident is None:
                    ident = float(tag[5:])
                elif tag.startswith('dv:f:') and ident is None:
                    ident = 1.0 - float(tag[5:])
                elif tag.startswith('ch:Z:'):
                    chain = tag[5:].split('.')[0]
            qlen[q] = int(t[1])
            recs.append((q, r, qs, qe, alen, ident, chain))
    return recs, qlen


def union_len(intervals):
    if not intervals:
        return 0
    intervals = sorted(intervals)
    total = 0
    cs, ce = intervals[0]
    for s, e in intervals[1:]:
        if s > ce:
            total += ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    return total + (ce - cs)


def main(path):
    recs, qlen = load(path)
    genome = sum(qlen.values())

    bp_per_pair = collections.defaultdict(collections.Counter)
    for q, r, qs, qe, a, i, c in recs:
        bp_per_pair[q][r] += a
    dom = {q: c.most_common(1)[0][0] for q, c in bp_per_pair.items()}

    all_cov = collections.defaultdict(list)
    diag_cov = collections.defaultdict(list)
    chain_refs = collections.defaultdict(set)
    id_weighted = 0.0
    id_bp = 0.0
    for q, r, qs, qe, a, i, c in recs:
        all_cov[q].append((qs, qe))
        if r == dom[q]:
            diag_cov[q].append((qs, qe))
        if i is not None:
            id_weighted += i * a
            id_bp += a
        if c is not None:
            chain_refs[(q, c)].add(r)

    total_bp = sum(r[4] for r in recs)
    diag_bp = sum(bp_per_pair[q][dom[q]] for q in bp_per_pair)
    diag_cov_sum = sum(union_len(v) for v in diag_cov.values())
    all_cov_sum = sum(union_len(v) for v in all_cov.values())

    stats = {
        'n': len(recs),
        'genome': genome,
        'total_bp': total_bp,
        'diag_bp': diag_bp,
        'precision': diag_bp / total_bp if total_bp else 0.0,
        'diag_cov': diag_cov_sum,
        'diag_cov_frac': diag_cov_sum / genome if genome else 0.0,
        'all_cov': all_cov_sum,
        'all_cov_frac': all_cov_sum / genome if genome else 0.0,
        'w_ident': id_weighted / id_bp if id_bp else 0.0,
        'cross_ref_chains': sum(1 for v in chain_refs.values() if len(v) > 1),
    }
    for k, v in stats.items():
        print(f"{k}={v:.4f}" if isinstance(v, float) else f"{k}={v}")


if __name__ == "__main__":
    main(sys.argv[1])
