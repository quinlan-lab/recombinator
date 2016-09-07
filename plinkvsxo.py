"""
We want to see if crossovers are associated with an increased rate of de novo
mutations. The problem is that not all crossovers will be assoicated with a
nearby de novo.

We iterate over each crossover in a family (mom+dad) and tabulate the number
of times we see a de novo within YY bases in either kid. This is the observed

To get the expected, iterate over each crossover in a family; we take all
families that do not have a crossover within XX bases and and tabulate the
number of times we see a denovo within YY bases.

We can then do chi-sq on obs, exp.

In the future, we can see how to optimize XX and YY.

"""
import glob
import random
import toolshed as ts
from interlap import InterLap
from collections import defaultdict
import scipy.stats as ss
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

extend = 0

fam_lookup = {d['sampleid']: d['sscfam'] for d in ts.reader('~/Data/ssc_519.ordered.ped')}

fams = list(set(fam_lookup.values()))
assert len(fams) == 519

xos = defaultdict(InterLap)

f = "results/2016_08_29-shapeit-2.5M/aggregated-filtered.bed"

for d in ts.reader(f):
    xos[d['chrom']].add((int(d['start']), int(d['end']), d['family_id']))
        #xos[d['chrom']].add((int(d['start']), int(d['end']), random.choice(fams)))


novos = defaultdict(InterLap)
header = "chrom start end _ _ _ _ mom dad kid _ _ _".split()

for d in ts.reader("/scratch/ucgd/lustre/u1006375/simons/519.plinkseq.denovos.bed",
        header=header):
    novos[d['chrom'][3:]].add((int(d['start']), int(d['end']), fam_lookup[d['kid']]))

seen = defaultdict(InterLap)
vals = [0, 0, 0, 0]


contingency = []

counts = defaultdict(int)

ors = []
pvs = []
for chrom in sorted(xos.keys()):

    for start, end, family_id in xos[chrom]:
        # find all overlaps in this region
        xo_hits = list(xos[chrom].find((start, end)))
        s, e = min(x[0] for x in xo_hits), max(x[1] for x in xo_hits)
        if (s, e) in seen[chrom]: continue
        seen[chrom].add((s, e))

        xo_hits = list(xos[chrom].find((s, e)))

        novo_hits = list(novos[chrom].find((s - extend, e + extend)))
        novo_fams = [x[2] for x in novo_hits]
        xo_fams = [x[2] for x in xo_hits]

        #xo_fams = set(xo_fams)
        #novo_fams = set(novo_fams)


        # number of kids with dn / total
        p_dn = len(novo_hits) / (519 * 2.0)
        # number of kids with xo / total
        p_xo = len(xo_hits) / (519 * 2.0)

        # number of kids with dn / total
        p_dn = len(novo_hits) / (519 * 1.0)
        # number of kids with xo / total
        p_xo = len(xo_hits) / (519 * 1.0)

        #
        p_dn_given_xo = sum(1.0 for f in novo_fams if f in xo_fams) / len(xo_fams)

        counts[(len(novo_fams), len(xo_fams), p_dn_given_xo > (p_dn * p_xo))] += 1

        if len(novo_hits) > 20 and len(xo_hits) > 20:
            #        | dn | not dn |
            # xo     |  ? |   ?    |
            # not xo |  ? |   ?    |
            dn_and_xo = sum(1 for f in novo_fams if f in xo_fams)
            not_dn_and_xo = sum(1 for f in xo_fams if not f in novo_fams)
            dn_and_not_xo = sum(1 for f in novo_fams if not f in xo_fams)
            neither = 2 * 519 - dn_and_xo - not_dn_and_xo - dn_and_not_xo

            odds, pv = ss.fisher_exact([
                [dn_and_xo, not_dn_and_xo],
                [dn_and_not_xo, neither]
                ])
            print("%.3f,%.4g\n%s\n%s\n" % (odds, pv, [dn_and_xo, not_dn_and_xo], [dn_and_not_xo, neither]))
            if odds < 1:
                odds = 1/-odds
            ors.append(odds)
            pvs.append(pv)

ys = []
sizes = []

fig, axes = plt.subplots(2)

xs = range(1, 51)
print "cutoff\ttrues\tfalses\tpval"
for cutoff in xs:

    trues = sum(c for dn, xo, c in counts if xo >= cutoff)
    false = sum(not c for dn, xo, c in counts if xo >= cutoff)
    pval = ss.binom_test(trues, false + trues)
    print "%d\t%d\t%d\t%.4g" % (cutoff, trues, false, pval)
    sizes.append(trues + false)

    ys.append(pval)

ys = -np.log10(ys)
sizes = np.array(sizes, dtype=float)
sizes /= sizes.mean()
sizes *= 30.0

axes[0].plot(xs, ys, marker=None, ls='--')
axes[0].scatter(xs, ys, marker='o', s=sizes)
axes[0].set_xlim(xs[0], xs[-1])
axes[0].set_ylim(ymin=0)
axes[0].set_xlabel("required number of crossovers")
axes[0].set_ylabel("-log10(p-value)")
axes[0].set_title("binomial p-value for P(DN|XO) > P(DN) * P(XO) varying number of XO events\npoints sized by number of denovos + number of xos at each cutoff")


ppvs = [p for o, p in zip(ors, pvs) if o > 1]
opvs = [p for o, p in zip(ors, pvs) if o <= 1]
ppvs = -np.log10(ppvs)
opvs = -np.log10(opvs)

axes[1].hist([ppvs, opvs], 20, label=['OR>1', 'OR<=1'], log='y')
axes[1].set_xlabel("-log10(fisher-p)")
axes[1].set_ylabel("count")
axes[1].set_ylim(ymin=0)
axes[1].set_xlim(xmin=0)
axes[1].legend()
plt.show()
