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
        T[d['chrom']].add((int(d['start']), int(d['end']), key_fmt.format(**d)))
    return T

def same(a, b):
    return a.startswith(b) or b.startswith(a)

def regionvs(n_samples, regions, query, region_fmt, query_fmt, extend=0):
    Q = mktree(query, query_fmt)
    R = mktree(regions, region_fmt)
    counts = defaultdict(int)

    for chrom in sorted(R):
        seen = InterLap()
        for start, end, region_key in R[chrom]:
            r_hits = list(R[chrom].find((start, end)))
            # make sure we don't look at the same region more than once.
            s, e = min(r[0] for r in r_hits), max(r[1] for r in r_hits)
            if (s, e) in seen: continue
            seen.add((s, e))

            # exapnd the selection to the full range:
            r_hits = list(R[chrom].find((s, e)))
            q_hits = list(Q[chrom].find((s - extend, e + extend)))

            pR = len(r_hits) / n_samples
            pQ = len(q_hits) / n_samples


            p_Q_given_R = sum(any(same(r[2], q[2]) for q in q_hits) for r in r_hits) / float(len(r_hits))

            # check if the joint probability is > than expected by chance.
            gt = p_Q_given_R > (pR * pQ)

            ret = {'pR': pR, 'nR': len(r_hits),
                   'pQ': pQ, 'nQ': len(q_hits),
                   'p_Q_given_R': p_Q_given_R,
                   'gt': gt,
                   'R_and_Q': sum(1.0 for q in q_hits if any(same(q[2], r[2]) for r in r_hits)),
                   'Q_and_not_R': sum(1.0 for q in q_hits if not any(same(q[2], r[2]) for r in r_hits)),
                   'R_and_not_Q': sum(1.0 for r in r_hits if not any(same(r[2], q[2]) for q in q_hits)),
                   }

            ret['neither_R_nor_Q'] = n_samples - sum((ret['R_and_not_Q'], ret['R_and_Q'], ret['Q_and_not_R']))
            yield ret


def run(args):
    n = args.n_samples
    r = ts.reader(args.regions)
    q = ts.reader(args.query)
    for d in regionvs(n, r, q, args.region_key, args.query_key):
        print(d)

if __name__ == "__main__":
    sys.exit(main())
