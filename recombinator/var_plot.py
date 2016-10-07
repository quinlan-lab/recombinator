from __future__ import absolute_import, print_function
import sys
import toolshed as ts
from os.path import basename
from matplotlib import pyplot as plt
import seaborn as sns
import math

def main(argv):
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--log", type=int, help="log the values.")
    p.add_argument("--column", required=True, help="column to plot")
    p.add_argument("--beds", nargs="+", help="bed file(s) from variant-info or denovo")
    args = p.parse_args(argv)
    run(args)

def tryfloat(n, log=False, vmax=100):
    try:
        v = float(n)
        if log:
            if v == 0: return vmax
            v = math.log(v, abs(log))
            if log < 0:
                v = -v
        return v
    except:
        raise
        return n

def run(args):
    if len(args.beds) == 1:
        just_one = True
    else:
        just_one = False
    for bed in args.beds:
        vals = []
        labels = []
        vals.append([tryfloat(d[args.column], args.log) for d in ts.reader(bed)])
        labels.append(basename(bed))
        if just_one:
            plt.hist(vals, label=labels, bins=20)
        elif not just_one:
            plt.hist(vals, label=labels, alpha=0.5, bins=20)
    
    plt.legend()
    plt.savefig("hist.png")

if __name__ == "__main__":
    if __package__ is None:
        from os import path
        sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
    import doctest
    doctest.testmod()
    main(sys.argv[1:])
