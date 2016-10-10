from __future__ import print_function, absolute_import
import sys
import time
from collections import OrderedDict

from peddy import Ped
from cyvcf2 import VCF, Writer

import numpy as np
import scipy.stats as ss
from . import denovo

def variant_ok(v, HET, exclude=None, min_mean_depth=20, max_mean_depth=150, min_pval=0.05,
        min_variant_qual=30, min_call_rate=0.99):
    """
    returns True if this is a good candidate for phasing.
    exclude should be the interval tree for the chromosome matching the
    variant.
    Recommended to use LCR:
    and sedgupds:
    http://humanparalogy.gs.washington.edu/build37/data/GRCh37GenomicSuperDup.tab
    with: cat GRCh37GenomicSuperDup.tab | awk '{a=NF-1; if(NR > 1 && $a < 0.02) { print substr($1, 4)"\t"$2"\t"$3  }}'
    (zcat data/LCR-hs37d5.bed.gz && cat segdup.ltp02.bed) | sort -k1,1 -k2,2n | bedtools merge -i stdin | bgzip -c > data/LCR-and-segdup.bed.gz
    for exclude.
    """
    if not denovo.variant_prefilter(v, min_variant_qual):
        return False
    if len(v.ALT) > 1:
        return False

    if v.call_rate < min_call_rate: continue

    m = np.mean(v.gt_depths)
    if not (min_mean_depth < m < max_mean_depth): return False

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
    p.add_argument("--denovos", help="if specified, de novos will be written to this file")
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

    fh_denovos = open(args.denovos, "w") if args.denovos else None
    if fh_denovos is not None:
        samples_lookup = {v: i for i, v in enumerate(vcf.samples)}
        kids = [k for k in ped.samples() if k.mom is not None and k.dad is not None]

    kept, HET = 0, vcf.HET
    t0 = time.time()
    n_dn = 0

    for i, v in enumerate(vcf(args.chrom) if args.chrom else vcf, start=1):
        if i % 100000 == 0:
            secs = time.time() - t0
            persec = i / float(secs)
            extra = (" called %d de-novos" % n_dn) if fh_denovos is not None else ""
            print("kept: %d of %d (%.3f%%). %.2f/sec.%s" % (kept, i,
                                                            100. * float(kept) / i,
                                                            persec, extra),
                  file=sys.stderr)
        if not denovo.variant_prefilter(v, 10):
            continue
        if fh_denovos is not None:
            d = denovo.get_denovo(v, samples_lookup, kids, max_alts_in_parents=1, exclude=exclude)
            if d is not None:
                denovo.write_denovo(d, fh_denovos)
                n_dn += 1

        if not variant_ok(v, HET, exclude): continue
        wtr.write_record(v)
        kept += 1
    secs = time.time() - t0
    persec = i / float(secs)
    extra = (" called %d de-novos" % n_dn) if fh_denovos is not None else ""
    print("kept: %d of %d (%.3f%%). %.2f/sec.%s" % (kept, i,
                                                    100. * float(kept) / i,
                                                    persec, extra),
          file=sys.stderr)
    wtr.close()


if __name__ == "__main__":
    main(sys.argv[1:])
    if __package__ is None:
        from os import path
        sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
