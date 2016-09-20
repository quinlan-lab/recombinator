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

def mktree(iterable, sample_key):
    T = defaultdict(InterLap)
    for d in iterable:
        #if 'parent_sex' in d and d['parent_sex'] != 'female': continue
        if d['chrom'] in ('X', 'Y'): continue
        T[d['chrom']].add([int(d['start']), int(d['end']), d[sample_key]])
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
    if randomize:
        R = copy.deepcopy(R)
        shuffle(R)

    for chrom in R:
        seen = InterLap()
        for start, end, sample_id in R[chrom]:
            if end - start > size_cutoff: continue

            q_hits = set(x[2] for x in Q[chrom].find((start - extend, end + extend)))
            n += int(sample_id in q_hits)
    return n


def _get_ovl(args):
    if args[2]:
        return [get_overlap_counts(*args) for _ in range(10)]
    else:
        return get_overlap_counts(*args)


def enrichment(regions, query, region_key, query_key, extend=10000,
             simulations=1000,
             size_cutoff=80000,
             figpath=None):
    R = mktree(regions, region_key)
    Q = mktree(query, query_key)

    obs = get_overlap_counts(Q, R)

    colors = sns.color_palette()

    res = []
    # divide simulations by 10 because we do 10 in _get_ovl to aid
    # parallelization.
    with ProcessPoolExecutor() as pool:
        for r in pool.map(_get_ovl, ((Q, R, True, size_cutoff, extend) for
            i in range(simulations/10)), chunksize=5):
            res.extend(r)
    print(obs)
    p = (1.0 + sum(r >= obs for r in res)) / float(len(res))
    print(p)

    fig, ax = plt.subplots(1)
    plt.title("Co-occurence")
    ax.hist(res, 25, label="expected")
    ax.axvline(x=obs, label="observed", color=colors[1], lw=3)
    ax.text(0.7, 0.92, "p: %.4g" % p, transform=ax.transAxes)
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
