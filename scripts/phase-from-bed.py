"""
an example script of how to phase de novo variants using GATK ReadBackedPhasing
"""

from __future__ import print_function
import toolshed as ts
import shutil
import os
import atexit
import glob

vcf = "/uufs/chpc.utah.edu/common/home/u6000771/Data/519FamiliesUnrecal_snp-recal_indel-recal.vcf.gz"
fasta = "/uufs/chpc.utah.edu/common/home/u6000771/bcbio/genomes/Hsapiens/g1k_v37_decoy/seq/g1k_v37_decoy.fa"
gatk = "/uufs/chpc.utah.edu/common/home/u6000771/bcbio/anaconda/opt/gatk-3.6/GenomeAnalysisTK.jar"
bam_glob = "/scratch/ucgd/lustre/u1006375/simons/de_novo/data/brent_minibams/{sample_id}.chr{chrom}.{start}*.sorted.bam"

def rm(path):
    try:
        os.unlink(path)
    except OSError:
        pass

window = 500
for d in ts.reader(1):

    samples = "{sample_id},{paternal_id},{maternal_id}".format(**d)
    region = "{chrom}:{start}-{end}".format(chrom=d['chrom'], start=int(d['start']) - window, end=int(d['end']) + window)

    # extract just the region and samples we need.
    ivcf = "{sample_id}.vcf".format(**d)
    cmd = "|bcftools view -a -c1 -m2 -M2 -s {samples} {vcf} {region} > {ivcf}".format(samples=samples, region=region, vcf=vcf, ivcf=ivcf)
    phasedvcf = "{sample_id}.phased.vcf".format(**d)

    atexit.register(rm, ivcf)
    atexit.register(rm, phasedvcf)
    atexit.register(rm, ivcf + ".idx")
    atexit.register(rm, phasedvcf + ".idx")

    bam = glob.glob(bam_glob.format(**d))
    if len(bam) != 1:
        print(bam_glob.format(**d))
        continue
    bam = bam[0]

    list(ts.nopen(cmd))
    #list(ts.nopen("""|bash rbp.sh {sample_id}".format(**d)))
    list(ts.nopen("|java -Xmx4G -jar {gatk} \
     -T ReadBackedPhasing -R {fasta} -I {bam} --variant {ivcf} -o {phasedvcf}".format(**locals())))

    found = False
    if not os.path.exists(phasedvcf):
        continue

    posns = []
    for line in ts.nopen(phasedvcf):
        if line[0] == "#": continue
        if "|" in line:
            found = True
            posns.append(int(line.split("\t", 4)[1]))
    if found:
        dst = "{sample_id}-{chrom}-{start}-{end}.phased.vcf".format(**d)
        shutil.move(phasedvcf, dst)
        region = "{chrom}:{start}-{end}".format(chrom=d['chrom'],
                start=min(posns)-10, end=max(posns)+10)
        print("samtools tview -p {region} {bam} {fasta}".format(**locals()))

    rm(phasedvcf + ".idx")
    rm(phasedvcf)
    rm(ivcf)
    rm(ivcf + ".idx")
