import argparse
import sys
import itertools as it
from collections import defaultdict
from toolshed import reader

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--bg", dest='bedgraph', required=True)
    args = p.parse_args()
    run(args)


def report(p, c):
    if c is not None:
        print '\t'.join(str(s) for s in [p[0], int(p[2])-1, int(c[1])+1])
    else:
        print '\t'.join(str(s) for s in [p[0], int(p[2])-1, int(p[2])])


def run(args):
    """
    1   19941   8199263 1
    1   8213867 64090377    2
    1   64096093    165191702   1
    1   165193661   249238596   2

    yields

    1   8199262 8213868
    1   64090376    64096094
    1   165191701   165193662
    """
    last_block = {}
    prev = None
    curr = None
    for l in open(args.bedgraph):
        curr = l.strip().split()
        if prev is None:
            pass
        else:
            if curr[0] != prev[0]:
                pass
            else:
                if curr[3] != prev[3]:
                    report(prev, curr)
        prev = curr



if __name__ == "__main__":
    main()
