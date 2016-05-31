from __future__ import print_function
import numpy as np
from collections import Counter
from hmmlearn import hmm
import sys

np.random.seed(42)


def fit(obs, n=20, transmat=np.array([[0.999, 0.001], [0.001, 0.999]])):
    """
    >>> obs = [0] * 20 + [0, 1, 0, 1, 0] + [1] * 20
    >>> fit(obs)
    array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

    >>> obs = [1] * 20 + [0, 1, 0, 1, 0] + [0] * 20
    >>> fit(obs)
    array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])


    >>> obs = [0] + [1] * 20 + [0, 1, 0, 1, 0] + [0] * 20
    >>> fit(obs)
    array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    """
    aobs = np.asarray(obs)
    if aobs.max() == 2:
        aobs -= 1
    print(Counter(aobs), file=sys.stderr)

    obs = aobs.reshape(1, -1).T
    model = hmm.GaussianHMM(n_components=2, n_iter=100,
            init_params='mcs',
            params="mcs")
    model.transmat_ = transmat
    model.fit(obs)
    res = model.predict(obs)

    if np.abs(res - aobs).mean() > 0.5:
        res = 1 - res
    return res


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
    if res.failed > 0:
        print(res, file=sys.stderr)
        sys.exit(1)
    main()
