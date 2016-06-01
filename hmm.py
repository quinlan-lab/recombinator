from __future__ import print_function
import numpy as np
from collections import Counter
from hmmlearn import hmm
import pomegranate as po
import sys

np.random.seed(42)

eps = 1e-7

def fit(obs, noise_pct=0.3, eps=eps, pseudocount=300):
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
    >>> obs = [0, 1, 0, 1, 0] + [1] * 400
    >>> res = fit(obs)
    >>> np.all(res == 1)
    True
    """
    obs = np.asarray(obs)
    if obs.max() == 2:
        obs -= 1

    if noise_pct > 1:
        noise_pct /= 100.


    model = po.HiddenMarkovModel("crossover")

    # allow, e.g. 6% of sites in a true state 1 to appear as state 0
    di0 = po.DiscreteDistribution({0: 1 - noise_pct, 1: noise_pct})
    di1 = po.DiscreteDistribution({1: 1 - noise_pct, 0: noise_pct})

    dnoise = po.DiscreteDistribution({0: 0.5, 1: 0.5})

    i0 = po.State(di0, name="0")
    i1 = po.State(di1, name="1")
    inoise = po.State(dnoise, name="2")

    model.add_states([i0, i1, inoise])

    m = np.mean(obs[:200])

    model.add_transition(model.start, i0, m - 0.001)
    model.add_transition(model.start, i1, 1 - m - 0.001)
    model.add_transition(model.start, inoise, 0.002)


    model.add_transition(i0, i0, 1 - 2.0 * eps, pseudocount=pseudocount)
    model.add_transition(i0, i1, eps)
    model.add_transition(i0, inoise, eps)

    model.add_transition(i1, i1, 1 - 2.0 * eps, pseudocount=pseudocount)
    model.add_transition(i1, i0, eps)
    model.add_transition(i1, inoise, eps)

    model.add_transition(inoise, inoise, 1 - 20 * eps)
    model.add_transition(inoise, i0, 10 * eps)
    model.add_transition(inoise, i1, 10 * eps)

    model.bake()
    _, path = model.viterbi(obs)
    #return np.array([x[0] for x in path[1:]])
    return np.array([int(x[1].name) for x in path[1:]])

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
