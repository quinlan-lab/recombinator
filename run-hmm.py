from slurmpy import Slurm
import os

creds = {"account": "quinlan-kp", "partition": "quinlan-kp"}

fams = set([x.split()[0] for i, x in
        enumerate(open("/uufs/chpc.utah.edu/common/home/u6000771/Data/ssc_519.ped"))
        if i > 0])

base = "/uufs/chpc.utah.edu/common/home/u1007787/Projects/recombination/results/"
pre = """
set +e
set +o pipefail
set -o nounset
PATH=/scratch/ucgd/lustre/u6000771/gem/tools/bin:/scratch/ucgd/lustre/u6000771/gem/data/anaconda/bin/:$PATH:~u6000771/bin
BASE={base}
""".format(base=base)

tmpl = """
IN=$BASE/split/chr{chrom}-fam{fam}-{parent}.bed.gz
if [[ ! -f $IN ]]; then
(zcat $BASE/{chrom}-0-1000000.bed.gz| head -1; (zgrep -hw {fam} $BASE/{chrom}-*.bed.gz | grep -w {parent})) \\
        | bgzip -c > $IN
fi

zcat $IN \\
        | python hmm.py $BASE/crossovers/chr{chrom}-fam{fam}-{parent}.bed  \\
        | bgzip -c > $BASE/hmm/chr{chrom}-fam{fam}-{parent}.hmm.bed.gz &

"""

plot_tmpl = """
python plotter.py $BASE/hmm/chr{chrom}-fam{fam}-{parent}.hmm.bed.gz \\
        $BASE/hmm/chr{chrom}-fam{fam}-{parent}.hmm.png &
"""

jobs = []
plot_jobs = []

cmd = "sbatch"

for fam in fams:
    for parent in ("mom", "dad"):
        for chrom in range(1, 23) + ["X", "Y"]:

            png = "{base}/hmm/chr{chrom}-fam{fam}-{parent}.hmm.png".format(**locals())
            if os.path.exists(png):
                continue

            jobs.append(tmpl.format(**locals()))
            plot_jobs.append(plot_tmpl.format(**locals()))

            if len(jobs) >= 44:

                s = Slurm("hmm-%s-%s-%s" % (fam, parent, chrom), creds)
                script = pre + "\n".join(jobs) + "\nwait\n" + "\n".join(plot_jobs) + "\nwait"
                s.run(script, name_addition="", _cmd=cmd)
                jobs = []
                plot_jobs = []


if len(jobs):
    s = Slurm("hmm-%s-%s-%s" % (fam, parent, chrom), creds)
    script = pre + "\n".join(jobs) + "\nwait\n" + "\n".join(plot_jobs) + "\nwait"
    s.run(script, name_addition="", _cmd=cmd)
    jobs = []
    plot_jobs = []
