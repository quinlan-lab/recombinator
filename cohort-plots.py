from __future__ import print_function
import sys
import math
import itertools as it
import array
from collections import defaultdict

import numpy as np
from scipy import stats
from matplotlib import ticker, pyplot as plt
import seaborn as sns
sns.set_style('whitegrid')

def maxend(fbed):
    emax = 0
    for toks in (x.rstrip().split("\t") for x in open(fbed)):
        try:
            emax = max(int(toks[2]), emax)
        except ValueError:
            continue
    return emax

def main(fbed, ped=None, prefix=None):
    e = maxend(fbed)
    d = {'1': np.zeros(e, dtype=np.uint16),
         '2': np.zeros(e, dtype=np.uint16)}

    lens, last_chrom = [], None

    sample_counts = defaultdict(lambda: defaultdict(int))
    sex = {x[1]: x[4] for x in (l.split("\t", 5) for l in open(ped))}

    for i, toks in enumerate(x.rstrip().split("\t") for x in open(fbed)):
        if i == 0:
            header = toks
            parent_id = toks.index('parent_id')
            continue
        pid = toks[parent_id]
        sample_counts[sex[pid]][pid] += 1

        if last_chrom != toks[0]:
            if last_chrom is not None:
                report(last_chrom, d, lens, prefix)
                d['1'][:] = 0
                d['2'][:] = 0
            lens = []
            last_chrom = toks[0]

        s, e = int(toks[1]), int(toks[2])
        lens.append(e - s)
        d[sex[pid]][s:e] += 1
    report(last_chrom, d, lens, prefix)
    plot_sample_counts(sample_counts, prefix)


def fmtr(x, p):
      v, suf = "%.2f" % (x / 1000000.), "M"
      # strip trailing 0's and "."
      while v[-1] in '.0' and '.' in v:
          v = v[:-1]
      return v + suf

def absfmtr(y, p):
    v = str(abs(y))
    while v[-1] in '.0' and '.' in v:
        v = v[:-1]
    return v

def plot_sample_counts(sample_counts, prefix):
    fig, ax = plt.subplots(1)
    ax.hist([np.array(sample_counts['1'].values()) / 2.0,
             np.array(sample_counts['2'].values()) / 2.0],
             20, label=['male', 'female'])
    ax.set_xlabel('Recombinations per meiosis')
    ax.set_ylabel('Count', rotation='vertical')
    plt.legend()
    plt.savefig('%s.recombinations-per-parent.png' % prefix)
    plt.close()


def report(chrom, d, lens, prefix, file=sys.stdout, zcutoff=2.58):

    fig, ax = plt.subplots(1, figsize=(12, 3))
    current_palette = sns.color_palette()

    for k, (sex, sex_label) in enumerate((('1', 'male'), ('2', 'female'), ('3', 'both'))):
        xs, ys, zs = array.array('I'), array.array('I'), array.array('f')
        if sex_label == 'both':
            # we output the total count for male+female
            # only print for the both and only plot for male, female separate
            arr = np.abs(d.pop('1')) + np.abs(d.pop('2'))
        else:
            arr = d[sex]
        zscore = stats.zscore(np.abs(arr))
        diff, = np.where(arr[:-1] != arr[1:])
        diff += 1
        for i, posn in enumerate(diff):
            if i == len(diff) - 1: break
            v = arr[posn]
            if v == 0: continue

            end = diff[i+1]
            vals = arr[posn:end]
            assert len(set(vals)) == 1, (vals, i, len(diff))
            z = zscore[posn]
            print("%s\t%d\t%d\t%d\t%s\t%.3f" % (chrom, posn, end, vals[0],
                  sex_label, z), file=file)

            xs.extend((posn, end))
            ys.extend((vals[0], vals[0]))
            zs.extend((z, z))

        posn = diff[-1]
        if arr[posn] != 0:
            xs.extend((posn, posn + 1))
            ys.extend((arr[posn], arr[posn]))
            zs.extend((zscore[posn], zscore[posn]))

            print("%s\t%d\t%d\t%d\t%s\t%.3f" % (chrom, posn, posn + 1,
                  arr[posn], sex_label, zscore[posn]), file=file)

        if sex_label == "both": break
        xs, ys, = np.asarray(xs, dtype=int), np.asarray(ys, dtype=int)
        zs = np.asarray(zs)

        line = ys[zs > zcutoff].min()
        if k == 0:
            ys, line = -ys, -line

        ax.plot(xs, ys, '-', label=sex_label, color=current_palette[k])
        ax.axhline(y=line, ls='--', color='0.4')

    ax.get_xaxis().set_major_formatter(ticker.FuncFormatter(fmtr))
    ax.get_yaxis().set_major_formatter(ticker.FuncFormatter(absfmtr))
    ax.text(0.015, 0.88, "chromosome: " + chrom, transform=ax.transAxes)
    ax.legend(loc='upper right')
    ax.set_xlabel('Genomic position')
    ax.set_ylabel('Samples with crossover', rotation='vertical')

    plt.tight_layout()
    plt.savefig("%s%s.png" % (prefix.rstrip("."), chrom if prefix.endswith("/") else ("." + chrom)))
    plt.close(fig)

if __name__ == "__main__":
    import sys
    f = sys.argv[1]
    ped = sys.argv[2]
    prefix = sys.argv[3]
    main(f, ped, prefix)
