from __future__ import print_function, absolute_import
import sys
import time
from collections import OrderedDict

from peddy import Ped
from cyvcf2 import VCF, Writer

import numpy as np
import scipy.stats as ss

def variant_prefilter(v, min_variant_qual):

    if len(v.REF) > 4: return False
    if len(v.ALT) > 2 or "*" in v.ALT: return False
    if len(v.ALT[0]) > 3: return False
    if v.FILTER is not None: return False

    if v.QUAL < min_variant_qual: return False
    return True

def get_denovo(v, samples, kids, max_alts_in_parents=1,
        min_depth=8,
        min_allele_balance_p=0.1,
        HET=1):
    """
    v: cyvcf2.Variant
    samples: dictionary of sample: index in vcf.
    kids: list of peddy.Samples() that have parents
    max_alts_in_parents: max number of alternate reads that can appear in the parents.
    """

    if v.num_het > 2: return None
    if v.num_hom_alt > 2: return None

    ret = []
    gts = v.gt_types
    ref_depths, alt_depths = None, None
    for kid in kids:
        ki = samples[kid.sample_id]
        if gts[ki] != HET: continue
        # check that parents are hom-ref
        mi, di = samples[kid.mom.sample_id], samples[kid.dad.sample_id]
        for pi in (mi, di):
            if gts[pi] != 0: continue


        if alt_depths is None:
            ref_depths = v.gt_ref_depths
            alt_depths = v.gt_alt_depths

        if alt_depths[[mi, di]].sum() > max_alts_in_parents: continue

        kid_alt = alt_depths[ki]
        kid_ref = ref_depths[ki]
        if kid_ref + kid_alt < min_depth: continue

        if ref_depths[di] + alt_depths[di] < min_depth: continue
        if ref_depths[mi] + alt_depths[mi] < min_depth: continue

        # if there are too many alts outside this kid. skip
        asum = alt_depths.sum() - kid_alt
        if asum > len(samples) / 10.:
            continue

        # now do a statistical test for the same thing:
        palt = ss.binom_test([asum, asum + ref_depths.sum() - kid_ref], p=0.0004,
                alternative="greater")
        if palt < min_allele_balance_p: continue

        pab = ss.binom_test([kid_ref, kid_alt])
        if pab < min_allele_balance_p: continue
        ret.append(OrderedDict((
            ("chrom", v.CHROM),
            ("start", v.start),
            ("end", v.end),
            ("sample_id", kid.sample_id),
            ("family_id", kid.family_id),
            ("ref", v.REF),
            ("alt", ",".join(v.ALT)),
            ("pab", pab),
            ("palt", palt),
            ("paternal_id", kid.paternal_id),
            ("maternal_id", kid.maternal_id),
            ("kid_ref_depth", kid_ref),
            ("kid_alt_depth", kid_alt),
            ("mom_ref_depth", ref_depths[mi]),
            ("mom_alt_depth", alt_depths[mi]),
            ("dad_ref_depth", ref_depths[di]),
            ("dad_alt_depth", alt_depths[di]),
            ("cohort_alt_depth", asum),
            )))

    # shouldn't have multiple samples with same de novo.
    if len(ret) == 1: return ret[0]

def variant_ok(v, HET, exclude=None, min_mean_depth=20, min_pval=0.05, min_variant_qual=30):
    """
    returns True if this is a good candidate for phasing.
    exclude should be the interval tree for the chromosome matching the
    variant.
    """
    if not variant_prefilter(v, min_variant_qual):
        return False
    if len(v.ALT) > 1:
        return False

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
    p.add_argument("--denovos", help="if specified, de novos will be written to this file")
    p.add_argument("ped")
    p.add_argument("vcf")

    run(p.parse_args(argv))

def write_denovo(d, fh_denovos):
    if fh_denovos.tell() == 0:
        print("#" + "\t".join(d.keys()), file=fh_denovos)
    print("\t".join(map(str, d.values())), file=fh_denovos)

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

    fh_denovos = open(args.denovos, "w") or None
    if fh_denovos is not None:
        samples_lookup = {v: i for i, v in enumerate(vcf.samples)}
        kids = [k for k in ped.samples() if k.mom is not None and k.dad is not None]

    HET = vcf.HET
    kept = 0

    t0 = time.time()

    for i, v in enumerate(vcf(args.chrom) if args.chrom else vcf, start=1):
        if not variant_prefilter(v, 10):
            continue
        if fh_denovos is not None:
            d = get_denovo(v, samples_lookup, kids, max_alts_in_parents=1)
            if d is not None:
                write_denovo(d, fh_denovos)

        if not variant_ok(v, HET, exclude): continue
        wtr.write_record(v)
        kept += 1
        if i % 100000 == 0:
            secs = time.time() - t0
            persec = i / float(secs)
            print("kept: %d of %d (%.3f%%). %.2f/sec." % (kept, i, 100. * float(kept) / i, persec),
                  file=sys.stderr)
    secs = time.time() - t0
    persec = i / float(secs)
    print("kept: %d of %d (%.3f%%). %.2f/sec." % (kept, i, 100. * float(kept) / i, persec),
          file=sys.stderr)
    wtr.close()


if __name__ == "__main__":
    main(sys.argv[1:])
    if __package__ is None:
        from os import path
        sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
