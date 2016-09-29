from __future__ import print_function, absolute_import
import sys
import time

from peddy import Ped
from cyvcf2 import VCF, Writer

import numpy as np
import scipy.stats as ss

def variant_ok(v, HET, exclude=None, min_mean_depth=20, min_pval=0.05, min_variant_qual=30):
    """
    returns True if this is a good candidate for phasing.
    exclude should be the interval tree for the chromosome matching the
    variant.
    """

    if len(v.REF) > 3: return False
    if len(v.ALT) > 1: return False
    if len(v.ALT[0]) > 3: return False
    if v.FILTER is not None: return False

    if v.QUAL < min_variant_qual: return False
    if np.mean(v.gt_depths) < min_mean_depth: return False

    hets, = np.where(v.gt_types == HET)
    if len(hets) == 0: return False


    alts, depths = v.gt_alt_depths, v.gt_depths

    vab = alts[hets] / depths[hets].astype(float)
    # test that all hets are in these bounds.
    if np.any(~((vab > 0.3) & (vab < 0.7))): return False

    if any(ss.binom_test(alts[k], depths[k]) < min_pval for k in hets):
        return False
    if exclude is not None and 0 != len(exclude[v.CHROM].search(v.start, v.end)):
        return False
    return True

def main(argv):
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument("--chrom", help="extract only this chromosome", default=None)
    p.add_argument("--exclude", help="regions to exclude (e.g. LCRs)")
    p.add_argument("--min-mean-depth", default=20, type=int)
    p.add_argument("ped")
    p.add_argument("vcf")

    run(p.parse_args(argv))

def run(args):

    vcf = VCF(args.vcf)
    ped = Ped(args.ped)
    psamples = set([x.sample_id for x in ped.samples()])
    samples = [x for x in vcf.samples if x in psamples]

    if args.exclude:
        from .recombinator import read_exclude
        exclude = read_exclude(args.exclude, args.chrom)
    else:
        exclude = None

    vcf = VCF(args.vcf, gts012=True, samples=samples)
    wtr = Writer("-", vcf)
    HET = vcf.HET
    kept = 0

    t0 = time.time()

    for i, v in enumerate(vcf(args.chrom) if args.chrom else vcf, start=1):
        if not variant_ok(v, HET, exclude): continue
        wtr.write_record(v)
        kept += 1
        if i % 100000 == 0:
            secs = time.time() - t0
            persec = i / float(secs)
            print("kept: %d of %d (%.3f%%). %.2f/sec." % (kept, i, 100. * float(kept) / i, persec),
                  file=sys.stderr)
    print("kept: %d of %d (%.3f%%). %.2f/sec." % (kept, i, 100. * float(kept) / i, persec),
          file=sys.stderr)
    wtr.close()



if __name__ == "__main__":
    main(sys.argv[1:])
    if __package__ is None:
        import sys
        from os import path
        sys.path.append( path.dirname( path.dirname( path.abspath(__file__) ) ) )
