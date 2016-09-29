"""
Evaluate co-ocurrence between `regions` and `query`
"""
from __future__ import print_function

import sys
from argparse import ArgumentParser
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import random
import copy

from interlap import InterLap
from peddy import Ped
import toolshed as ts

try:
    from itertools import izip
    zip = izip
    range = xrange
except ImportError:
    pass

from matplotlib import pyplot as plt
import seaborn as sns
sns.set_style('whitegrid')

def main(args=sys.argv[1:]):

    p = ArgumentParser(__doc__)
    p.add_argument("--region-key", default="sample_id", help="column-header from regions to use for sample")
    p.add_argument("--query-key", default="sample_id", help="column-header from query to user for sample")
    p.add_argument("--figure", default="xodn.png", help="path indicating where to save the figure")
    p.add_argument("--simulations", type=int, default=1000, help="number of shufflings to compare to observed")
    p.add_argument("--extend", type=int, default=1, help="extend regions on either size by this amount before checking for overlaps.")
    p.add_argument("--size-cutoff", type=int, default=80000, help="exclude regions greater than this length")

    p.add_argument("regions", metavar="REGIONS")
    p.add_argument("query", metavar="QUERY")
    a = p.parse_args(args)
    return run(a)

def mktree(iterable, sample_key, size_cutoff):
    T = defaultdict(InterLap)
    for d in iterable:
        #if 'parent_sex' in d and d['parent_sex'] != 'female': continue
        if d['chrom'] in ('X', 'Y'): continue
        s, e = int(d['start']), int(d['end'])
        if e - s > size_cutoff: continue
        T[d['chrom']].add([s, e, d[sample_key]])
    return T

def get_overlap_counts(Q, R, size_cutoff=80000, extend=10000, do_print=False):
    n = 0
    a_pairs, b_pairs = [], []
    for chrom in R:
        for start, end, sample_id in R[chrom]:
            if end - start > size_cutoff: continue

            q_hits = list(Q[chrom].find((start - extend, end + extend)))
            q_samples = [x[2] for x in q_hits]
            a_pairs.extend([sample_id] * len(q_samples))
            b_pairs.extend(q_samples)
            if do_print:
                for st, en, smp in q_hits:
                    if smp != sample_id: continue
                    print("%s\t%d\t%d\t%s\t%d\t%d" % (chrom, start, end, sample_id, st, end))

            n += int(sample_id in set(q_samples))
    return n, a_pairs, b_pairs


def enrichment(regions, query, region_key, query_key, extend=10000,
               simulations=1000,
               size_cutoff=180000,
               figpath=None):
    R = mktree(regions, region_key, size_cutoff)
    Q = mktree(query, query_key, size_cutoff)

    obs, aa_pairs, bb_pairs = get_overlap_counts(Q, R)
    print("observed:", obs)
    print("n-pairs:", len(aa_pairs))
    print(sum(a == b for a, b in zip(aa_pairs, bb_pairs)))
    res = []
    import time
    import numpy as np

    lookup = {p: i for i, p in enumerate(set(aa_pairs + bb_pairs))}

    # faster to compair ints than strings so we convert here.
    a_pairs = np.array([lookup[p] for p in aa_pairs], dtype=np.int)
    b_pairs = np.array([lookup[p] for p in bb_pairs], dtype=np.int)
    obs = (a_pairs == b_pairs).sum()

    print("obs2:", obs)
    t0, t1 = time.time(), time.time()
    for i in range(simulations):
        if i > 0 and i % 10000 == 0:
            print(i, time.time() - t1, res[-1])
            t1 = time.time()
        np.random.shuffle(a_pairs)
        res.append((a_pairs == b_pairs).sum())
    ngt = sum(r >= obs for r in res)

    p2p5, p50, p97p5 = np.percentile(res, [5, 50, 95])

    p = (1.0 + ngt) / float(1 + len(res))
    print(ngt, p)
    print(time.time() - t0)

    colors = sns.color_palette()

    enriched = obs / float(p50)
    e2p5, e97p5 = obs / float(p2p5), obs / float(p97p5)

    fig, ax = plt.subplots(1)
    plt.title("Co-occurrence")
    ax.hist(res, min(max(res), 25), label="expected")
    ax.axvline(x=obs, label="observed", color=colors[1], lw=3)
    #ax.text(0.66, 0.92, "p: %.3g (1 + %d) / (1 + %d)" % (p, ngt, len(res)), transform=ax.transAxes)
    ax.text(0.60, 0.92, "p: %.3g\nFC:%.2f (%.2f-%.2f)"
            % (p, enriched, e97p5, e2p5), transform=ax.transAxes)
    ax.set_xlabel("Number of overlaps")
    ax.set_ylabel("Count")
    plt.legend(loc='upper left')
    plt.savefig(figpath)

def run(args):
    r = ts.reader(args.regions)
    q = ts.reader(args.query)
    enrichment(r, q, args.region_key, args.query_key,
            simulations=args.simulations,
            extend=args.extend,
            size_cutoff=args.size_cutoff,
            figpath=args.figure)

if __name__ == "__main__":
    sys.exit(main())
