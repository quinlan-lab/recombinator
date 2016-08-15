from slurmpy import Slurm
from collections import defaultdict
import os

creds = {"account": "quinlan-kp", "partition": "quinlan-kp"}
PED = "/uufs/chpc.utah.edu/common/home/u6000771/Data/ssc_519.ped"
DATE = "2016_07_21"
DATE = "2016_08_12-shapeit-sites"
DATE = "2016_08_12-shapeit-filtered"

fams = set([x.split()[0] for i, x in enumerate(open(PED)) if i > 0])

fam_kids = defaultdict(list)
keys = "family_id sample_id dad_id mom_id".split()
for toks in (l.rstrip().split() for l in open(PED) if l[0] != '#'):
    d = dict(zip(keys, toks))
    if d['dad_id'] in ('0', '-9', ''): continue
    if d['mom_id'] in ('0', '-9', ''): continue
    fam_kids[d['family_id']].append(d['sample_id'])

base = "/scratch/ucgd/lustre/u6000771/Projects/src/recombinator/results/{DATE}/".format(**locals())

pre = """
set +e
set +o pipefail
set -o nounset
PATH=/scratch/ucgd/lustre/u6000771/gem/tools/bin:/scratch/ucgd/lustre/u6000771/gem/data/anaconda/bin/:$PATH:~u6000771/bin
BASE={base}

mkdir -p $BASE/hmm/
mkdir -p $BASE/crossovers/
""".format(base=base)

tmpl = """

# writes .bed and .filtered.bed for crossovers.
zcat $BASE/recomb/{chrom}/fam{fam}/{chrom}-*.{parent}-{kid}.bed.gz \\
        | python hmm.py $BASE/crossovers/chr{chrom}-fam{fam}-{parent}-{kid}  \\
        | bgzip -c > $BASE/hmm/chr{chrom}-fam{fam}-{parent}-{kid}.hmm.bed.gz &

"""

plot_tmpl = """
python plotter.py $BASE/hmm/chr{chrom}-fam{fam}-{parent}-{kid}.hmm.bed.gz \\
        $BASE/crossovers/chr{chrom}-fam{fam}-{parent}-{kid}.filtered.bed \\
        $BASE/hmm/chr{chrom}-fam{fam}-{parent}-{kid}.hmm.png &
"""

jobs = []
plot_jobs = []

cmd = "sbatch"

for fam in fams:
    kids = fam_kids[fam]
    if kids == []:
        sys.stderr.write("WARNING! no kids found for family %s\n" % fam)
        continue
    for kid in kids:
        for parent in ("mom", "dad"):
            #for chrom in range(1, 23) + ["X", "Y"]:
            for chrom in [22]:

                png = "{base}/hmm/chr{chrom}-fam{fam}-{parent}-{kid}.hmm.png".format(**locals())
                if os.path.exists(png):
                    print(png)
                    continue

                jobs.append(tmpl.format(**locals()))
                plot_jobs.append(plot_tmpl.format(**locals()))

                if len(jobs) >= 44:

                    s = Slurm("hmm-%s-%s-%s-%s" % (fam, parent, kid, chrom), creds)
                    script = pre + "\n".join(jobs) + "\nwait\n" + "\n".join(plot_jobs) + "\nwait"
                    s.run(script, name_addition="", _cmd=cmd)
                    jobs = []
                    plot_jobs = []


if len(jobs):
    s = Slurm("hmm-%s-%s-%s-%s" % (fam, parent, kid, chrom), creds)
    script = pre + "\n".join(jobs) + "\nwait\n" + "\n".join(plot_jobs) + "\nwait"
    s.run(script, name_addition="", _cmd=cmd)
    jobs = []
    plot_jobs = []
