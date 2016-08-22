from __future__ import print_function
import sys
import os
import itertools as it
from operator import itemgetter
from collections import Counter

import numpy as np
import toolshed as ts
import pomegranate as po

np.random.seed(42)

eps = 1e-12
MIN_INFORMATIVE_SITES = 200

def fit(obs, noise_pct=0.08, eps=eps, pseudocount=0):
    """
    >>> obs = [0, 1, 0, 1, 0] + [1] * 40
    >>> res = fit(obs)
    >>> np.all(res == 1)
    True
    """
    obs = np.asarray(obs)
    if obs.max() == 2:
        obs -= 1
    #print("obs:", list(obs), file=sys.stderr)

    if noise_pct > 1:
        noise_pct /= 100.

    model = po.HiddenMarkovModel("crossover")

    # allow, e.g. 6% of sites in a true state 1 to appear as state 0
    di0 = po.DiscreteDistribution({0: 1 - noise_pct, 1: noise_pct})
    di1 = po.DiscreteDistribution({1: 1 - noise_pct, 0: noise_pct})

    dnoise = po.DiscreteDistribution({0: 0.5, 1: 0.5})

    i0 = po.State(di0, name="0")
    i1 = po.State(di1, name="1")
    #  inoise = po.State(dnoise, name="2")

    model.add_states([i0, i1])

    m = np.mean(obs[:MIN_INFORMATIVE_SITES])
    #print("mean:", m, file=sys.stderr)

    model.add_transition(model.start, i0, 1 - m)
    model.add_transition(model.start, i1, m)
    #model.add_transition(model.start, inoise, 0.002)


    model.add_transition(i0, i0, 1 - 1.0 * eps, pseudocount=pseudocount)
    model.add_transition(i0, i1, eps)
    #model.add_transition(i0, inoise, eps)

    model.add_transition(i1, i1, 1 - 1.0 * eps, pseudocount=pseudocount)
    model.add_transition(i1, i0, eps)
    #model.add_transition(i1, inoise, eps)

    #model.add_transition(inoise, inoise, 1 - 20 * eps)
    #model.add_transition(inoise, i0, 10 * eps)
    #model.add_transition(inoise, i1, 10 * eps)

    model.bake()
    _, path = model.viterbi(obs)
    #return np.array([x[0] for x in path[1:]])
    return np.array([int(x[1].name) for x in path[1:]])

def main(args=None):

    import toolshed as ts
    if args is None:
        if len(sys.argv) == 1:
            print("usage: %s <x-over-bed> <input>\nif input is not specified, stdin is assumed" % sys.argv[0])

        args = sys.argv[1:]
    if len(args) == 1:
        args.append("-")

    rows = []
    xfh = ts.nopen(args[0] + ".bed", mode="w")
    xfh.write("chrom\tstart\tend\tparent\tfamily_id\tstate-change\tinformative-sites\tmean_val\n")
    for d in ts.reader(args[1], header='ordered'):
        if float(d.get('call_rate', 1)) < 0.95: continue
        if d['start'] == 'start': continue
        v = int(d['same'])
        rows.append((d['chrom'], int(d['start']), v, d))
    rows.sort()

    if len(rows) == 0:
        return
    vals = fit([r[2] for r in rows])
    #print(vals, file=sys.stderr)
    #1/0

    last_start, last_state = None, None
    prev_start = None
    n = 0

    state_vals = []

    for i, row in enumerate((r[-1] for r in rows)):
        if i == 0:
            print("#" + "\t".join(row.keys() + ['hmm-state']))
            parent_col, kid_col = [k for k in row.keys() if k[-1] == '*']
        print("%s\t%d" % ("\t".join(row.values()), vals[i]))
        if last_state is None:
            last_state = vals[i]
            last_start = row['start']
        n += 1

        if vals[i] != last_state:
            d = row.copy()
            d['start'] = prev_start
            d['end'] = row['start']
            d['parent'] = d[parent_col]
            d['state-change'] = "%d-%d" % (last_state, vals[i])
            d['informative-sites'] = str(n)
            d['mean_val'] = sum(state_vals) / float(len(state_vals))
            xfh.write("{chrom}\t{start}\t{end}\t{parent}\t{family_id}\t{state-change}\t{informative-sites}\t{mean_val}\n".format(**d))
            last_start = row['start']
            last_state = vals[i]
            state_vals = []
            n = 0

        state_vals.append(int(row['same']))
        prev_start = row['start']
    xfh.close()
    xo_filter(xfh.name, xfh.name.rsplit(".", 1)[0] + ".filtered.bed")


def xo_filter(iname, oname, delta=0.05,
        min_informative_sites=MIN_INFORMATIVE_SITES):
    """
    filter x-overs to minimum informative sites and purity (delta).
    merges info for x-os that were previously separated by a spurios x-o
    """

    last = None
    vals = [x for x in list(ts.reader(iname, header="ordered")) if
            int(x['informative-sites']) > min_informative_sites and (float(x['mean_val']) < delta or float(x['mean_val']) > 1 - delta)]
    if len(vals) == 0:
        if os.path.exists(oname):
            os.unlink(oname)
        return

    ofh = ts.nopen(oname, "w")
    ofh.write("\t".join(vals[0].keys()) + "\n")

    for _, grp in it.groupby(vals, itemgetter('state-change')):
        grp = list(grp)
        grp[-1]['informative-sites'] = str(sum(int(g['informative-sites']) for g in grp))
        ofh.write("\t".join(grp[-1].values()) + "\n")
    ofh.close()


if __name__ == "__main__":
    """
    import doctest
    res = doctest.testmod()
    print(res, file=sys.stderr)
    if res.failed > 0:
        sys.exit(1)
    """
    main()
