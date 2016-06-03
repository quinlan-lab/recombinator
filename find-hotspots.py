"""
given the small set of actual crossovers called by hmm.py,
find hotspots from the files on stdin in bin_sizes given by sys.argv[1].
"""
from collections import defaultdict
import sys


counts = defaultdict(lambda: defaultdict(float))

window_size = int(sys.argv[1])


header = None
for i, line in enumerate(sys.stdin):
    if i == 0:
        header = line.rstrip().split("\t")
        continue
    d = dict(zip(header, line.rstrip().split()))
    if d['start'] == 'start': continue  # header
    if '2' in d['state-change']: continue  # unknown state


    try:
        chrom = int(d['chrom'])
    except:
        chrom = d['chrom']

    sidx, _ = divmod(int(d['start']), window_size)
    eidx, _ = divmod(int(d['end']), window_size)

    counts[chrom][sidx] += 0.5
    counts[chrom][eidx] += 0.5


for chrom in sorted(counts):
    dchrom = counts[chrom]
    for idx in sorted(dchrom):
        count = dchrom[idx]
        start  = idx * window_size
        end = (idx + 1) * window_size
        print("{chrom}\t{start}\t{end}\t{count}".format(**locals()))
