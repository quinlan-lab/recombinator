"""
convert a VCF and a PED to a .dat, .map, and a .ped for use by merlin.

"""
import sys
import gzip
import os
import re

from cyvcf2 import VCF
import numpy as np
from peddy import Ped
import tempfile
import atexit

# assert all(a == b for a, b in zip(samples, vcf.samples)), (samples[:10], vcf.samples[:10])

def main(args=sys.argv[1:]):
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--region", help="optional region, e.g. '1:112335-117798' or just '11'", default="")
    p.add_argument("--min-maf", help="minimum maf to consider", default=0.03, type=float)
    p.add_argument("ped")
    p.add_argument("vcf")
    p.add_argument("output_prefix")

    a = p.parse_args(args)
    run(a)

def run(args):

    p = Ped(args.ped)
    samples = [s.sample_id for s in p.samples()]
    vcf = VCF(args.vcf, samples=samples, gts012=True)

    #fhs = [open(tempfile.mkstemp()[1], "w") for _ in samples]
    fhs = [[] for _ in samples]
    #for fh in fhs:
    #    atexit.register(os.unlink, fh.name)

    dat = open("%s.dat" % args.output_prefix, "w")
    dat.write("A	affection\n")

    fmap = open("%s.map" % args.output_prefix, "w")
    fmap.write("CHROMOSOME\tMARKER\tPOSITION\n")

    # yes, '0 0' is missing
    #l = np.array(["1/ 1\n", "1/ 2\n", "2/ 2\n", "0/ 0\n"])
    l = np.array(["1/ 1", "1/ 2", "2/ 2", "0/ 0"])
    #l = np.array(["0/ 0\n", "0/ 1\n", "1/ 1\n", "x x\n"])
    #l = np.array(["A A\n", "A C\n", "C C\n", "x x\n"])

    k = 0
    for v in vcf(args.region):
        if v.QUAL < 20: continue
        if v.call_rate < 0.9: continue
        if v.aaf < args.min_maf: continue
        if np.mean(v.gt_depths) < 20: continue
        if np.mean(v.gt_depths) > 400: continue
        if np.mean(v.gt_quals) < 20: continue

        if not v.FILTER in ("PASS", None): continue
        if len(v.ALT) > 1: continue

        k += 1
        if k == 1600000:
            break
        marker = "%s_%d" % (v.CHROM, v.start + 1)
        if k % 10000 == 0:
            sys.stderr.write("%d @ %s:%d\n" % (k, v.CHROM, v.start + 1))
        dat.write("M	%s\n" % marker)

        bases = l[v.gt_types]
        bases[v.gt_quals < 20] = l[-1]
        fmap.write("%s\t%s\t%.8f\n" % (v.CHROM, marker, v.start / 1000000.0))
        for i, fh in enumerate(fhs):
            fh.append(bases[i])
    dat.close()
    fmap.close()


    #for fh in fhs: fh.close()

    merlin_ped = open("%s.ped" % args.output_prefix, "w")

    lsample = {s.sample_id: s for s in p.samples()}

    for i, sample_id in enumerate(vcf.samples):
        sample = lsample[sample_id]
        s = str(sample.affected and 2 or 1)
        # need replace twice because python won't do overlapping replacements.
        s = ("\t".join(str(sample).split()[:5]) + "\t" + s).replace("	-9	", "\t0\t").replace("	-9	", "\t0\t")
        #lines = [x.strip() for x in open(fhs[i].name).read().split("\n")]
        lines = fhs[i]
        merlin_ped.write("%s\t%s\n" % (s, "\t".join(lines)))
    merlin_ped.close()

if __name__ == "__main__":
    main()
