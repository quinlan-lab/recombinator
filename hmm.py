from __future__ import print_function
import numpy as np
from collections import Counter
from hmmlearn import hmm
import pomegranate as po
import sys

np.random.seed(42)

eps = 1e-8

def fit(obs, noise_pct=0.2):
    """
    #>>> obs = [0] * 200 + [0, 1, 0, 1, 0] + [1] * 200
    #>>> fit(obs)
    array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

    #>>> obs = [1] * 20 + [0, 1, 0, 1, 0] + [0] * 20
    #>>> fit(obs)
    array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0,
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    #>>> obs = [0, 1, 0, 1, 0] + [1] * 400
    #>>> res = fit(obs)
    #>>> np.all(res == 1)
    #True
    """
    obs = np.asarray(obs)
    if obs.max() == 2:
        obs -= 1

    if noise_pct > 1:
        noise_pct /= 100.


    model = po.HiddenMarkovModel("crossover")

    # allow, e.g. 6% of sites in a true state 1 to appear as state 0
    di0 = po.DiscreteDistribution({1: 1 - noise_pct, 0: noise_pct})
    di1 = po.DiscreteDistribution({0: 1 - noise_pct, 1: noise_pct})

    i0 = po.State(di0, name="0")
    i1 = po.State(di1, name="1")

    model.add_states(i0, i1)

    m = np.mean(obs[:50])

    model.add_transition(model.start, i0, m)
    model.add_transition(model.start, i1, 1 - m)

    model.add_transition(i0, i0, 1 - eps)
    model.add_transition(i1, i1, 1 - eps)
    model.add_transition(i0, i1, eps)
    model.add_transition(i1, i0, eps)

    model.bake()
    _, path = model.viterbi(obs)
    #model.fit(obs)
    v = np.array([x[0] for x in path[1:]])
    if np.abs(v - obs).mean() > 0.5:
        v = 1 - v
    return v

def main(args=None):

    import toolshed as ts
    if args is None:
        if len(sys.argv) > 1:
            args = sys.argv[1:]
        else:
            args = ["-"]

    rows = []
    for d in ts.reader(args[0], header='ordered'):
        if d['start'] == 'start': continue
        v = int(d['same(1)_diff(2)'])
        rows.append((d['chrom'], int(d['start']), v, d))
    rows.sort()

    vals = fit([r[2] for r in rows])

    for i, row in enumerate((r[-1] for r in rows)):
        if i == 0:
            print("\t".join(row.keys() + ['hmm-state']))
        print("%s\t%d" % ("\t".join(row.values()), vals[i]))



if __name__ == "__main__":
    import doctest
    res = doctest.testmod()
    print(res, file=sys.stderr)
    if res.failed > 0:
        sys.exit(1)

    main()
