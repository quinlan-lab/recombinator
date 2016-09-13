# python regionvs.py --region-key "{kid_id}||{parent_id}" --query-key
# "{sample_id}||" results/2016_08_29-shapeit-2.5M/aggregated-filtered.bed
# plinkseq.denovos.bed 1038
"""
Evaluate overlap betwee `regions` and `query`
"""
from __future__ import print_function

import sys
from argparse import ArgumentParser
from collections import defaultdict
import random

from interlap import InterLap
from peddy import Ped
import toolshed as ts

import numpy as np
import seaborn as sns
sns.set_style('whitegrid')
from matplotlib import pyplot as plt
import scipy.stats as ss

def main(args=sys.argv[1:]):

    p = ArgumentParser(__doc__)
    p.add_argument("--region-key", required=True, help="python format string using headers from regions file, e.g {family_id},{sample_id}")
    p.add_argument("--query-key", required=True, help="python format string using headers from query file, e.g {family_id},{sample_id}. share prefix with regions")

    p.add_argument("regions", metavar="REGIONS")
    p.add_argument("query", metavar="QUERY")
    p.add_argument("n_samples", type=float, metavar="N")

    a = p.parse_args(args)
    return run(a)

def mktree(iterable, key_fmt):
    T = defaultdict(InterLap)
    for d in iterable:
        #if 'parent_sex' in d and d['parent_sex'] != 'female': continue
        if d['chrom'] in ('X', 'Y'): continue
        T[d['chrom']].add([int(d['start']), int(d['end']), key_fmt.format(**d)])
    return T

def same(a, b):
    return a == b
    return a.startswith(b) or b.startswith(a)

def shuffle(itree, z=[0]):
    vals = []
    for k in itree:
        # shuffle vals across all chroms.
        vals.extend(v[2] for v in itree[k]._iset)
    """
    import collections
    with open("t.%d.txt" % z[0], 'w') as fh:
        for k, v in collections.Counter(vals).items():
            fh.write("%s\t%d\n" % (k, v))
        z[0] += 1

    """
    random.shuffle(vals)
    rset = list(set(vals))
    start = 0
    for k in itree:
        tree = itree[k]
        assert isinstance(tree, InterLap)
        #vals = list(set([v[2] for v in tree._iset]))
        for i, v in enumerate(tree._iset):
            #tree._iset[i][2] = vals[i + start]
            tree._iset[i][2] = random.choice(rset)
        start += len(tree._iset)

def regionvs(n_samples, regions, query, region_fmt, query_fmt, extend=500,
             randomize=False):
    R = mktree(regions, region_fmt)
    Q = mktree(query, query_fmt)
    if randomize:
        shuffle(R)
        shuffle(Q)

    counts = defaultdict(int)
    n_samples = float(n_samples)

    for chrom in sorted(R):
        seen = InterLap()
        for start, end, region_key in R[chrom]:
            if end - start > 100000: continue
            """
            r_hits = list(R[chrom].find((start, end)))
            # make sure we don't look at the same region more than once.
            s, e = min(r[0] for r in r_hits), max(r[1] for r in r_hits)
            if (s, e) in seen: continue
            seen.add((s, e))
            """
            s, e = start, end
            if (s, e) in seen: continue
            seen.add((s, e))

            # expand the selection to the full range:
            r_hits = set([x[2] for x in R[chrom].find((s, e)) if x[1] - x[0] < 100000])
            q_hits = set([x[2] for x in Q[chrom].find((s - extend, e + extend)) if x[1] - x[0] < 100000])

            pR = len(r_hits) / n_samples
            pQ = len(q_hits) / n_samples

            if len(q_hits) > 0:
                if isinstance(q_hits, set):
                    qr_same = float(len(q_hits.intersection(r_hits)))
                else:
                    qr_same = sum(any(same(r, q) for r in r_hits) for q in q_hits)
                p_Q_given_R = qr_same / float(len(r_hits))
            else:
                qr_same = 0.0
                p_Q_given_R = 0.0

            # check if the joint probability is > than expected by chance.
            gt = p_Q_given_R > pQ

            ret = {'pR': pR, 'nR': len(r_hits),
                   'pQ': pQ, 'nQ': len(q_hits),
                   'p_Q_given_R': p_Q_given_R,
                   'region': '%s:%d-%d' % (chrom, s, e),
                   'gt': gt,
                   'R_and_Q': qr_same,
                   'Q_and_not_R': sum(1.0 for q in q_hits if not any(same(q[2], r[2]) for r in r_hits)),
                   'R_and_not_Q': sum(1.0 for r in r_hits if not any(same(r[2], q[2]) for q in q_hits)),
                   }

            ret['neither_R_nor_Q'] = n_samples - sum((ret['R_and_not_Q'], ret['R_and_Q'], ret['Q_and_not_R']))
            yield ret


def run(args):
    n = args.n_samples
    r = ts.reader(args.regions)
    q = ts.reader(args.query)
    res = list(regionvs(n, r, q, args.region_key, args.query_key))

    _xs, ys, sz, xs = range(1, 51), [], [], []
    diffs = []
    for cutoff in _xs:
        trues = sum(d['gt'] for d in res if d['nR'] >= cutoff)
        false = sum(not d['gt'] for d in res if d['nR'] >= cutoff)
        if trues + false == 0: break
        xs.append(cutoff)
        pval = ss.binom_test([trues, false], alternative='greater')
        print("cutoff:%d" % cutoff, "trues:", trues, 'falses:', false, 'p:', pval)
        ys.append(pval)
        sz.append(trues + false)
    diffs = [d['p_Q_given_R'] - d['pQ'] for d in res if d['nQ'] > 2]

    ys = -np.log10(ys)

    sz = np.array(sz, dtype=float)
    sz /= sz.mean()
    fig, axs = plt.subplots(3, figsize=(10, 10))
    ax = axs[0]
    ax.scatter(xs, ys, s=sz * 50)
    ax.plot(xs, ys, ls='--', marker=None)
    ax.set_xlim(xmin=xs[0], xmax=xs[-1])
    ax.set_ylim(ymin=0)
    ax.set_ylabel("-log10(p)")
    ax.set_xlabel("cutoff")
    #sns.stripplot(x='x', y='log2enrichment', ax=axs[1], data=df, jitter=0.1, color="0.3")


    ax = axs[1]
    ax.hist(diffs, 50)
    ax.set_xlabel("p(Q|R) - p(Q)")

    ax = axs[2]
    ax.scatter(
        #[d['p_Q_given_R'] for d in res],
        #[d['pQ'] * d['pR'] for d in res],
        [d['pQ'] for d in res],
        [d['pR'] for d in res],
        s=[d['nR'] for d in res])
    #ax.set_xlabel('p(Q|R)')
    #ax.set_ylabel('p(Q) * p(R)')
    ax.set_xlabel('p(Q)')
    ax.set_ylabel('p(R)')

    plt.show()
    plt.savefig("/uufs/chpc.utah.edu/common/home/u6000771/public_html/xxn.png")

if __name__ == "__main__":
    sys.exit(main())
