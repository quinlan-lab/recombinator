from __future__ import print_function
import sys
import math
import itertools as it
import array
from collections import defaultdict
import toolshed as ts

import numpy as np
from scipy import stats
import warnings

warnings.simplefilter('ignore')
from matplotlib import ticker, pyplot as plt
import seaborn as sns
sns.set_style('whitegrid')

cytoband_colors = dict(
    gpos100=(0,0,0),
    gpos=(0,0,0),
    gpos75=(100,100,100),
    gvar=(220,220,220),
    gpos66=(160,160,160),
    gpos50=(200,200,200),
    gpos33=(210,210,210),
    gpos25=(200,200,200),
    gneg=(255,255,255),
    acen=(217,47,39),
    stalk=(100,127,164),
)

for k, v in cytoband_colors.items():
    cytoband_colors[k] = (v[0] / 255., v[1] / 255., v[2] / 255.)

def maxend(fbed):
    emax = 0
    for toks in (x.rstrip().split("\t") for x in open(fbed)):
        try:
            emax = max(int(toks[2]), emax)
        except ValueError:
            continue
        except IndexError:
            print("bad line: %s", toks)
            continue
    return emax

def read_cytobands(cytobands):
    icytobands = defaultdict(list)
    if cytobands is None: return defaultdict(lambda: None)
    for toks in (l.rstrip().split() for l in ts.nopen(cytobands)):
        print(toks)
        chrom = toks[0].replace('chr', '')
        icytobands[chrom].append((chrom, int(toks[1]), int(toks[2]), cytoband_colors[toks[4]]))
    return icytobands

def run(xobeds, ped=None, prefix=None, min_size=0, max_size=sys.maxint,
        blocks=False, cytobands=None):
    sample_counts = defaultdict(lambda: defaultdict(int))
    sex = {x[1]: x[4] for x in (l.split("\t", 5) for l in open(ped))}
    cytobands = read_cytobands(cytobands)

    fhrates = open("%s.%s-rates.txt" % (prefix, "gc" if blocks else "xo"), "w")

    for xobed in sorted(xobeds):
        print(xobed)
        lens, last_chrom = [], None
        e = maxend(xobed)
        d = {'1': np.zeros(e, dtype=np.uint16),
             '2': np.zeros(e, dtype=np.uint16)}

        xfh = open(xobed)
        header = next(xfh).rstrip().split("\t")
        parent_id = header.index('parent_id')
        left_block = header.index('left-block')

        for i, toks in enumerate(sorted(x.rstrip().split("\t") for x in xfh)):

            if blocks:
                se = toks[left_block].split(":")[1].split("-")
                s, e = int(se[0]), int(se[1])
            else:
                s, e = int(toks[1]), int(toks[2])
            if not (min_size < e - s < max_size):
                print("skipping due to size:", e - s)
                continue

            pid = toks[parent_id]
            sample_counts[sex[pid]][pid] += 1

            if last_chrom != toks[0]:
                if last_chrom is not None:
                    print("switching chrom: %s->%s" % (last_chrom, toks[0]),
                          file=sys.stderr)
                    report_xo(last_chrom, d, lens, prefix, file=fhrates, cyto=cytobands[last_chrom])
                    d = {'1': np.zeros(e, dtype=np.uint16),
                         '2': np.zeros(e, dtype=np.uint16)}
                lens = []
                last_chrom = toks[0]

            lens.append(e - s)
            d[sex[pid]][s:e] += 1

        # get the last right block.
        if blocks and last_chrom is not None:
            se = toks[header.index('right-block')].split(":")[1].split("-")
            s, e = int(se[0]), int(se[1])
            if (min_size < e - s < max_size):
                sample_counts[sex[pid]][pid] += 1
                pid = toks[parent_id]
                lens.append(e - s)
                d[sex[pid]][s:e] += 1

        report_xo(last_chrom, d, lens, prefix, file=fhrates, cyto=cytobands[last_chrom])
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


def report_xo(chrom, d, lens, prefix, file=sys.stdout, zcutoff=2.58, cyto=None):

    fig, axes = plt.subplots(3, 1,
                gridspec_kw={'wspace': 0, 'hspace': 0, 'height_ratios': [15,
                    0.001 if cyto is None else 1 , 15]},
                sharex=True)
    ax = axes[0]
    plt.subplots_adjust(wspace=0, hspace=0)

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
            if posn > 0 and arr[posn-1] == 0:
                xs.append(posn-1)
                ys.append(0)
                zs.append(0)

            end = diff[i+1]
            vals = arr[posn:end]
            assert len(set(vals)) == 1, (vals, i, len(diff))
            z = zscore[posn]
            print("%s\t%d\t%d\t%d\t%s\t%.3f" % (chrom, posn, end, vals[0],
                  sex_label, z), file=file)

            xs.extend((posn, end))
            ys.extend((vals[0], vals[0]))
            zs.extend((z, z))
            if end+1 < len(arr) and arr[end+1] == 0:
                xs.append(end + 1)
                ys.append(0)
                zs.append(0)

        if len(diff) == 0:
            continue

        posn = diff[-1]
        if arr[posn] != 0:
            xs.extend((posn, posn + 1))
            ys.extend((arr[posn], arr[posn]))
            zs.extend((zscore[posn], zscore[posn]))

            print("%s\t%d\t%d\t%d\t%s\t%.3f" % (chrom, posn, posn + 1,
                  arr[posn], sex_label, zscore[posn]), file=file)

        if sex_label == "both": break
        xs, ys = np.asarray(xs, dtype=int), np.asarray(ys, dtype=int)
        zs = np.asarray(zs)

        line = ys[zs > zcutoff].min()
        if k == 0:
            ys, line = -ys, -line
            ax = axes[2]
        else:
            ax = axes[0]

        ax.plot(xs, ys, '-', label=sex_label, color=current_palette[k])
        ax.axhline(y=line, ls='--', color='0.4')

    ax = axes[0]
    ax.get_xaxis().set_major_formatter(ticker.FuncFormatter(fmtr))
    ax.get_yaxis().set_major_formatter(ticker.FuncFormatter(absfmtr))
    ax.text(0.015, 0.88, "chromosome: " + chrom, transform=ax.transAxes)
    ax.legend(loc='upper right')
    axes[1].set_ylabel('Samples with crossover', rotation='vertical')
    axes[1].yaxis.set_label_position("right")

    axes[2].set_xlabel('Genomic position')
    plt.draw()
    vals = ["%.0f" % abs(x.get_position()[1]) for x in axes[2].get_yticklabels()]
    axes[2].set_yticklabels(vals)
    axes[2].legend(loc='lower right')

    if cyto:
        ax = axes[1]
        for chrom, start, end, color in cyto:
            ax.axvspan(start, end, facecolor=color, edgecolor=color, lw=0, alpha=0.6)
        ax.grid(b=False, which='both')
        sns.despine(ax=axes[1], left=True, right=True)
    else:
        sns.despine(ax=axes[0], bottom=True, left=False, right=False, top=False)

    axes[1].set_xticks([])
    axes[1].set_yticks([])
    sns.despine(ax=axes[1], top=True, bottom=True)

    plt.tight_layout()
    plt.savefig("%s%s.png" % (prefix.rstrip("."), chrom if prefix.endswith("/") else ("." + chrom)))
    plt.close(fig)

def main(args):
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--ped", required=True, help="ped file")
    p.add_argument("--prefix", required=True, help="where to put images")
    p.add_argument("--blocks", action='store_true', default=False)
    p.add_argument("--min-size", default=0, type=int, help="mininum crossover size in bases to use.")
    p.add_argument("--max-size", default=10000, type=int, help="maximum crossover size in bases to use. if --blocks is specified, this si the maximum block-size to use.")
    p.add_argument("--cytobands", help="optional file of cytobands to plot")
    p.add_argument("crossover_files", nargs="+", help=".crossover.bed files")
    a = p.parse_args(args)
    run(a.crossover_files, a.ped, a.prefix, min_size=a.min_size,
        max_size=a.max_size, blocks=a.blocks, cytobands=a.cytobands)

if __name__ == "__main__":
    main(sys.argv[1:])
