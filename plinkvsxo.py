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
from scipy.stats import hypergeom
import scipy.stats as ss
import numpy as np

extend = 0

fam_lookup = {d['sampleid']: d['sscfam'] for d in ts.reader('~/Data/ssc_519.ordered.ped')}

fams = list(set(fam_lookup.values()))
assert len(fams) == 519

xos = defaultdict(InterLap)

for f in glob.glob("results/2016_08_29-shapeit-2.5M/crossovers/*.filtered.bed"):
    for d in ts.reader(f):
        xos[d['chrom']].add((int(d['start']), int(d['end']), d['family_id']))
        #xos[d['chrom']].add((int(d['start']), int(d['end']), random.choice(fams)))


novos = defaultdict(InterLap)
header = "chrom start end _ _ _ _ mom dad kid _ _ _".split()

for d in ts.reader("/scratch/ucgd/lustre/u1006375/simons/519.plinkseq.denovos.bed",
        header=header):
    novos[d['chrom'][3:]].add((int(d['start']), int(d['end']), fam_lookup[d['kid']]))
    #novos[d['chrom'][3:]].add((int(d['start']), int(d['end']), random.choice(fams)))

seen = defaultdict(InterLap)
# now we have the denovos and the xos in interval trees, we count the correspondence.
vals = [0, 0, 0, 0]

contingency = []

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


        #novo_fams = [random.choice(fams) for _ in novo_fams]
        #xo_fams = [random.choice(fams) for _ in xo_fams]

        xo_fams = set(xo_fams)
        nof = set(novo_fams)
        if len(nof) < 20: continue

        no_xo_fams = set(fams) - xo_fams

        #obs_shared = sum(1.0 for f in xo_fams if f in nof)
        #exp_shared = sum(1.0 for f in xo_fams if not f in nof)

        xo_and_de_novo = sum(1.0 for f in xo_fams if f in nof)
        no_xo_and_de_novo = sum(1.0 for f in no_xo_fams if f in nof)

        contingency.append((xo_and_de_novo, no_xo_and_de_novo, len(xo_fams), len(no_xo_fams)))
        print(contingency[-1])

