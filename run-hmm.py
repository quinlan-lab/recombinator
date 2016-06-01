from slurmpy import Slurm

creds = {"account": "quinlan-kp", "partition": "quinlan-kp"}

fams = set([x.split()[0] for i, x in
        enumerate(open("/uufs/chpc.utah.edu/common/home/u6000771/Data/ssc_519.ped"))
        if i > 0])

pre = """
set -eo pipefail -o nounset
PATH=/scratch/ucgd/lustre/u6000771/gem/tools/bin:/scratch/ucgd/lustre/u6000771/gem/data/anaconda/bin/:$PATH:~u6000771/bin
"""

tmpl = """
(zcat results/{chrom}-0-1000000.bed.gz| head -1; zgrep -hw {fam} results/{chrom}-*.bed.gz) \\
        | grep {parent} \\
        | python hmm.py  \\
        | bgzip -c > results/hmm/chr{chrom}-fam{fam}-{parent}.hmm.bed.gz &
"""

plot_tmpl = """
python plotter.py results/hmm/chr{chrom}-fam{fam}-{parent}.hmm.bed.gz \\
        results/hmm/chr{chrom}-fam{fam}-{parent}.hmm.png &
"""

jobs = []
plot_jobs = []

cmd = "sbatch"

for fam in fams:
    for parent in ("mom", "dad"):
        for chrom in range(1, 23) + ["X", "Y"]:

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
