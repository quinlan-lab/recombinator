import matplotlib
import sys
if len(sys.argv) > 2:
    matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt
import toolshed as ts

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
last_o = int(row['same(1)_diff(2)'])
last_h = int(row['hmm-state'])

xs, ys = [], []
for row in rdr:

    s, o, h = int(row['start']), int(row['same(1)_diff(2)']), int(row['hmm-state'])
    ys.append(o/2.0-0.5)
    xs.append(s)

    if o != last_o:
        axes[1].axvspan(xmin=last_ostart, xmax=s, ymin=(last_o / 2.0) - 0.5, ymax=last_o / 2.0)
        last_ostart = s

    if h != last_h:
        print("h", last_h, last_hstart, s)
        ym=(last_h + 1) / 3.0
        axes[2].axvspan(xmin=last_hstart, xmax=s, ymin=ym-0.33, ymax=ym)
        last_hstart = s

    last_o = o
    last_h = h

axes[0].plot(xs, ys, 'k.')
axes[0].set_ylabel('state')
axes[0].set_yticks([0, 0.5])
axes[0].set_yticklabels(['state0', 'state1'])

ym = (last_h + 1) / 3.0
axes[1].axvspan(xmin=last_ostart, xmax=s, ymin=(last_o / 2.0) - 0.5, ymax=last_o / 2.0)
axes[2].axvspan(xmin=last_hstart, xmax=s, ymin=ym-0.33, ymax=ym)

axes[1].set_xlim(0, s)
axes[1].set_ylabel('state')
axes[2].set_ylabel('state')

axes[1].set_yticks([0.25, 0.75])
axes[1].set_yticklabels(['state0', 'state1'])

axes[2].set_yticks([0.33 - 0.12, 0.66 - 0.12, 0.99 - 0.12])
axes[2].set_yticklabels(['state0', 'state1', 'unknown'])
if len(sys.argv) > 2:
    plt.savefig(sys.argv[2])
else:
    plt.show()
