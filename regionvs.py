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
        if 'parent_sex' in d and d['parent_sex'] != 'female': continue
        T[d['chrom']].add((int(d['start']), int(d['end']), key_fmt.format(**d)))
    return T

def same(a, b):
    return a.startswith(b) or b.startswith(a)

def regionvs(n_samples, regions, query, region_fmt, query_fmt, extend=0):
    Q = mktree(query, query_fmt)
    R = mktree(regions, region_fmt)
    counts = defaultdict(int)
    n_samples = float(n_samples)

    for chrom in sorted(R):
        seen = InterLap()
        for start, end, region_key in R[chrom]:
            r_hits = list(R[chrom].find((start, end)))
            # make sure we don't look at the same region more than once.
            s, e = min(r[0] for r in r_hits), max(r[1] for r in r_hits)
            if (s, e) in seen: continue
            seen.add((s, e))

            # expand the selection to the full range:
            r_hits = list(R[chrom].find((s, e)))
            q_hits = list(Q[chrom].find((s - extend, e + extend)))

            pR = len(r_hits) / n_samples
            pQ = len(q_hits) / n_samples

            qr_same = float(sum(any(same(r[2], q[2]) for q in q_hits) for r in r_hits))
            p_Q_given_R = qr_same / float(len(r_hits))

            # check if the joint probability is > than expected by chance.
            gt = p_Q_given_R >= (pR * pQ)

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

    xs, ys, sz = range(1, 101), [], []
    for cutoff in xs:
        trues = sum(d['gt'] for d in res if d['nR'] >= cutoff)
        false = sum(not d['gt'] for d in res if d['nR'] >= cutoff)

        pval = ss.binom_test(trues, false + trues, alternative='greater')
        print(trues, false, pval)
        ys.append(pval)
        sz.append(trues + false)

    ys = -np.log10(ys)

    sz = np.array(sz, dtype=float)
    sz /= sz.mean()
    plt.scatter(xs, ys, s=sz * 50)
    plt.plot(xs, ys, ls='--', marker=None)
    plt.xlim(xmin=xs[0], xmax=xs[-1])
    plt.ylim(ymin=0)
    plt.ylabel("-log10(p)")
    plt.xlabel("cutoff")
    plt.show()
    plt.savefig("/uufs/chpc.utah.edu/common/home/u6000771/public_html/x.png")

if __name__ == "__main__":
    sys.exit(main())
