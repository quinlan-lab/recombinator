import argparse
import sys
import itertools as it
from collections import defaultdict
from toolshed import reader

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--min-size", dest='min_size', type=int, default=50000)
    p.add_argument("--bg", dest='bedgraph', required=True)
    args = p.parse_args()
    run(args)


def report(block):
    print '\t'.join(str(s) for s in [block['chrom'], block['start'], block['end'], block['state']])


def run(args):
    """
    chrom   start   end     same(1)_diff(2) dad_SSC00369    mom_SSC00470    template_SSC00468   sib_SSC00471
    1       14932   14933   1               G/A             G/G             G/G                 G/G
    1       15446   15447   1               A/G             A/A             A/A                 A/A
    1       16256   16257   1               G/C             G/G             G/C                 G/C
    1       19917   19918   2               G/C             G/G             G/C                 G/G
    1       19918   19919   2               G/A             G/G             G/A                 G/G
    1       19941   19942   1               G/C             G/G             G/C                 G/C
    1       20155   20156   1               C/T             C/C             C/T                 C/T
    1       20157   20158   1               A/C             A/A             A/A                 A/A
    1       20165   20166   1               A/G             A/A             A/A                 A/A
    """
    last_block = {}
    prev = None
    curr = None
    for l in reader(args.bedgraph):
        curr = l
        if prev is None:
            last_block['chrom'] = curr['chrom']
            last_block['start'] = int(curr['start'])
            last_block['state'] = curr['same(1)_diff(2)']
            prev = curr

        else:
            if curr['chrom'] != prev['chrom']:
                report(last_block)
            else:
                if curr['same(1)_diff(2)'] != prev['same(1)_diff(2)']:
                    if last_block['end'] - last_block['start'] > args.min_size:
                        report(last_block)
                    last_block = {}
                    last_block['start'] = int(curr['start'])
                    last_block['chrom'] = curr['chrom']
                    last_block['state'] = curr['same(1)_diff(2)']
                    last_block['end'] = int(curr['end'])
                else:
                    last_block['end'] = int(curr['end'])
        prev = curr

    # cleanup of last block
    report(last_block)




if __name__ == "__main__":
    main()
