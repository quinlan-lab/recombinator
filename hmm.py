from __future__ import print_function
import numpy as np
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


if __name__ == "__main__":
    import doctest
    res = doctest.testmod()
    if res.failed > 0:
        print(res, file=sys.stderr)
        sys.exit(1)
