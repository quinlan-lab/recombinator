import matplotlib
import os
import sys
matplotlib.use('Agg')
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
import toolshed as ts


# http://stackoverflow.com/questions/13728392/moving-average-or-running-mean
def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0)).astype(float)
    return (cumsum[N:] - cumsum[:-N]) / N

fig, axes = plt.subplots(3, sharex=True, figsize=(12, 4))
axes[0].set_title("raw")
axes[1].set_title("x-over")
axes[2].set_title("hmm")

name = sys.argv[1].split("/")[-1].split(".")[0]
#fig.set_title(name)
fig.suptitle(name)

try:
    rdr = ts.reader(1)
    row = next(rdr)
except StopIteration:
    sys.stderr.write("%s was empty\n" % sys.argv[1])
    sys.exit(0)
last_hstart = last_ostart = int(row['start'])
last_o = int(row['same'])
last_h = int(row['hmm-state'])

xs, ys = [], []
for row in rdr:

    s, o, h = int(row['start']), int(row['same']), int(row['hmm-state'])
    ys.append(o)
    xs.append(s)

    if o != last_o:
        axes[1].axvspan(xmin=last_ostart, xmax=s, ymin=last_o / 2.0, ymax=last_o
                / 2.0 + 0.5, ec='none')
        last_ostart = s

    if h != last_h:
        ym=(last_h + 1) / 2.0
        axes[2].axvspan(xmin=last_hstart, xmax=s, ymin=ym-0.5, ymax=ym, ec='none')
        last_hstart = s

    last_o = o
    last_h = h

axes[0].plot(xs, ys, 'k.')


N = 40
ym = running_mean(ys, N)
# so we can see the line.
ym[ym == 1] = 0.98
ym[ym == 0] = 0.02
print(len(xs), len(ym))

axes[0].plot(xs[N-1:], ym, 'k--')



axes[0].set_ylabel('state')
axes[0].set_yticks([0, 0.5])
axes[0].set_yticklabels(['state0', 'state1'])

ym = (last_h + 1) / 2.0
axes[1].axvspan(xmin=last_ostart, xmax=s, ymin=last_o / 2.0, ymax=last_o / 2.0 +
        0.5, ec='none')
axes[2].axvspan(xmin=last_hstart, xmax=s, ymin=ym-0.5, ymax=ym, ec='none')

axes[1].set_xlim(0, s)
axes[1].set_ylabel('state')
axes[2].set_ylabel('state')

axes[1].set_yticks([0.25, 0.75])
axes[1].set_yticklabels(['state0', 'state1'])

axes[2].set_yticks([0.25, 0.75])
axes[2].set_yticklabels(['state0', 'state1'])

if len(sys.argv) > 3:
    # we got a list of filtered x-overs
    if os.path.exists(sys.argv[2]):
        for d in ts.reader(sys.argv[2]):
            for ax in axes:
                ax.axvspan(xmin=int(d['start']), xmax=int(d['end']), color='#ff3333', ymax=1.1)

    plt.savefig(sys.argv[3])
else:
    plt.savefig(sys.argv[2])
