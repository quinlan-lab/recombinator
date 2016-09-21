"""
Evaluate co-ocurrent between `regions` and `query`
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
    p.add_argument("--extend", type=int, default=10000, help="extend regions on either size by this amount before checking for overlaps.")
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

def shuffle(itree, z=[0]):
    vals = []
    for k in itree:
        # shuffle vals across all chroms.
        vals.extend(v[2] for v in itree[k]._iset)
    random.shuffle(vals)
    #rset = list(set(vals))
    start = 0
    for k in itree:
        tree = itree[k]
        assert isinstance(tree, InterLap)
        #vals = list(set([v[2] for v in tree._iset]))
        for i, v in enumerate(tree._iset):
            tree._iset[i][2] = vals[i + start]
            #tree._iset[i][2] = random.choice(rset)
        start += len(tree._iset)

def get_overlap_counts(Q, R, randomize=False, size_cutoff=80000, extend=10000):
    n = 0
    a_pairs, b_pairs = [], []
    if randomize:
        R = copy.deepcopy(R)
        shuffle(R)

    for chrom in R:
        seen = InterLap()
        for start, end, sample_id in R[chrom]:
            if end - start > size_cutoff: continue

            q_hits = [x[2] for x in Q[chrom].find((start - extend, end + extend))]
            a_pairs.extend([sample_id] * len(q_hits))
            b_pairs.extend(q_hits)

            n += int(sample_id in set(q_hits))
    return n, a_pairs, b_pairs


def _get_ovl(args):
    if args[2]:
        return [get_overlap_counts(*args) for _ in range(10)]
    else:
        return get_overlap_counts(*args)


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
    p = (1.0 + ngt) / float(1 + len(res))
    print(ngt, p)
    print(time.time() - t0)

    colors = sns.color_palette()

    """
    res = []
    # divide simulations by 10 because we do 10 in _get_ovl to aid
    # parallelization.
    with ProcessPoolExecutor() as pool:
        for r in pool.map(_get_ovl, ((Q, R, True, size_cutoff, extend) for
            i in range(simulations/10)), chunksize=5):
            res.extend(r)
    ngt = sum(r >= obs for r in res)
    p = (1.0 + ngt) / float(1 + len(res))
    print(p)
    """

    fig, ax = plt.subplots(1)
    plt.title("Co-occurrence")
    ax.hist(res, 25, label="expected")
    ax.axvline(x=obs, label="observed", color=colors[1], lw=3)
    #ax.text(0.66, 0.92, "p: %.3g (1 + %d) / (1 + %d)" % (p, ngt, len(res)), transform=ax.transAxes)
    ax.text(0.60, 0.92, "p: %.3g" % (p, ), transform=ax.transAxes)
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
