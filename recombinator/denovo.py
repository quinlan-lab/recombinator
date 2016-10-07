from __future__ import print_function, absolute_import, division
import sys
import time
from collections import OrderedDict

from peddy import Ped
from cyvcf2 import VCF, Writer

import numpy as np
import scipy.stats as ss

def tranche99(filt, cutoff=99.6):
    """
    return True if the tranche is below 99.6
    VQSRTrancheINDEL90.00to99.00
    """
    if filt is None: return True
    if filt[:4] != "VQSR": return False
    try:
        return float(filt.split("to")[1]) < cutoff
    except:
        return False


def variant_prefilter(v, min_variant_qual):

    if len(v.REF) > 4: return False
    if len(v.ALT) > 2 or "*" in v.ALT: return False
    if len(v.ALT[0]) > 6: return False
    if v.FILTER is not None and not tranche99(v.FILTER) : return False

    if v.QUAL < min_variant_qual: return False
    return True

def get_denovo(v, samples, kids, max_alts_in_parents=1,
        min_depth=5,
        max_mean_depth=400,
        min_allele_balance_p=0.05,
        min_qual=1,
        exclude=None,
        HET=1):
    """
    v: cyvcf2.Variant
    samples: dictionary of sample: index in vcf.
    kids: list of peddy.Samples() that have parents
    max_alts_in_parents: max number of alternate reads that can appear in the parents.
    min_depth: require all members of trio to have at least this depth.
    min_allele_balance_p: if the p-value for ref, alt counts is less than this, exclude
    min_depth_percentile: require all members of a trio to be in this percentile of depth.
    exclude: interval tree of regions to exclude.
    """

    if v.num_het > 2: return None
    if v.num_hom_alt > 1: return None

    ret = []
    gts = v.gt_types

    ref_depths, alt_depths = None, None
    for kid in kids:
        ki = samples[kid.sample_id]
        if gts[ki] != HET: continue
        # check that parents are hom-ref
        mi, di = samples[kid.mom.sample_id], samples[kid.dad.sample_id]
        for pi in (mi, di):
            if not (gts[pi] == 0 or gts[pi] == 3): continue

        if alt_depths is None:
            depths = v.format('AD', int)
            ref_depths = depths[:, 0]
            for k in range(1, depths.shape[1]):
                alt_depths = depths[:, k]
            #ref_depths = v.gt_ref_depths
            #alt_depths = v.gt_alt_depths
            alt_depths[alt_depths < 0] = 0
            ref_depths[ref_depths < 0] = 0

        if alt_depths[[mi, di]].sum() > max_alts_in_parents: continue

        kid_alt = alt_depths[ki]
        kid_ref = ref_depths[ki]

        # depth filtering.
        if kid_alt + kid_ref < min_depth + 1: continue
        if kid_alt < min_depth / 2: continue
        if ref_depths[di] + alt_depths[di] < min_depth - 1: continue
        if ref_depths[mi] + alt_depths[mi] < min_depth - 1: continue

        if np.mean(alt_depths + ref_depths) > max_mean_depth:
            continue

        # if there are too many alts outside this kid. skip
        alt_sum = alt_depths.sum() - kid_alt


        # this check is less stringent that the p-value checks below but
        # avoids some compute.
        if alt_sum > len(samples) * 0.01:
            continue

        # balance evidence in kid and alts in entire cohort.
        if alt_sum >= kid_alt: continue

        if alt_sum > 0:
            if kid_alt < 6: continue

        # via Tom Sasani.
        palt = ss.binom_test([alt_sum, ref_depths.sum() - kid_ref], p=0.0002,
                alternative="greater")
        if palt < min_allele_balance_p: continue

        pab = ss.binom_test([kid_ref, kid_alt])
        if pab < min_allele_balance_p: continue

        quals = v.gt_quals
        # TODO: check why some quals are 0 and if filtering on this improve
        # accuracy.
        #quals[quals < 0] == 0
        #if quals[ki] < 1 or quals[mi] < 1 or quals[di] < 1: continue

        # stricter settings with FILTER
        if v.FILTER is not None:
            # fewer than 1 alt per 500-hundred samples.
            if alt_sum > 0.002 * len(samples): continue
            # no alts in either parent.
            if alt_depths[[mi, di]].sum() > 0: continue

        ret.extend(variant_info(v, kid, samples, pab, palt))

    if exclude is not None and 0 != len(exclude[v.CHROM].search(v.start, v.end)):
        return None

    # shouldn't have multiple samples with same de novo.
    if len(ret) == 1: return ret[0]

def variant_info(v, kid, samples, pab=None, palt=None):

    quals = v.gt_quals
    ki, mi, di = samples[kid.sample_id], samples[kid.mom.sample_id], samples[kid.dad.sample_id]
    depths = v.format('AD', int)
    ref_depths = depths[:, 0]
    all_alts = depths[:, 1:]
    for k in range(all_alts.shape[1]):
        alt_depths = all_alts[:, k]

    #ref_depths, alt_depths, quals = v.gt_ref_depths, v.gt_alt_depths, v.gt_quals
        kid_ref, kid_alt = ref_depths[ki], alt_depths[ki]
        alt_sum = alt_depths.sum() - kid_alt
        if pab is None:
            pab = ss.binom_test([kid_ref, kid_alt])
        if palt is None:
            palt = ss.binom_test([alt_sum, ref_depths.sum() - kid_ref], p=0.0002,
                    alternative="greater")

        yield OrderedDict((
            ("chrom", v.CHROM),
            ("start", v.start),
            ("end", v.end),
            ("sample_id", kid.sample_id),
            ("family_id", kid.family_id),
            ("ref", v.REF),
            ("alt", v.ALT[k]),
            ("alt_i", "%d/%d" % (k, len(v.ALT))),
            ("filter", v.FILTER or "PASS"),
            ("pab", "%.3g" % pab),
            ("palt", "%.3g" % palt),
            ("paternal_id", kid.paternal_id),
            ("maternal_id", kid.maternal_id),
            ("kid_ref_depth", kid_ref),
            ("kid_alt_depth", kid_alt),
            ("kid_total_depth", kid_ref + kid_alt),
            ("mom_ref_depth", ref_depths[mi]),
            ("mom_alt_depth", alt_depths[mi]),
            ("mom_total_depth", ref_depths[mi] + alt_depths[mi]),
            ("dad_ref_depth", ref_depths[di]),
            ("dad_alt_depth", alt_depths[di]),
            ("dad_total_depth", ref_depths[di] + alt_depths[di]),
            ("kid_qual", quals[ki]),
            ("mom_qual", quals[mi]),
            ("dad_qual", quals[di]),
            #("p%d_depth" % min_depth_percentile, p5),
            ("depth_mean", "%.1f" % np.mean(ref_depths + alt_depths)),
            ("depth_std", "%.1f" % np.std(ref_depths + alt_depths)),
            ("qual_mean", "%.1f" % np.mean(quals)),
            ("call_rate", "%.3f" % v.call_rate),
            ("cohort_alt_depth", alt_sum),
            ))


def denovo(v, samples, kids, max_alts_in_parents=1,
        min_depth=5,
        max_mean_depth=400,
        min_allele_balance_p=0.05,
        min_depth_percentile=3,
        min_qual=1,
        exclude=None,
        HET=1):
    if not variant_prefilter(v, 10): return None
    return get_denovo(v, samples, kids, max_alts_in_parents=max_alts_in_parents,
            min_depth=min_depth,
            max_mean_depth=max_mean_depth,
            min_allele_balance_p=min_allele_balance_p,
            min_qual=min_qual,
            exclude=exclude,
            HET=HET)


def write_denovo(d, fh, _line={}):
    if _line.get(fh.name, 0) == 0:
        _line[fh.name] = 1
        print("#" + "\t".join(d.keys()), file=fh)
    print("\t".join(map(str, d.values())), file=fh)

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

    samples_lookup = {v: i for i, v in enumerate(vcf.samples)}
    kids = [k for k in ped.samples() if k.mom is not None and k.dad is not None]

    n_dn = 0
    for i, v in enumerate(vcf(args.chrom) if args.chrom else vcf, start=1):
        if i % 100000 == 0:
            print(" called %d de-novos out of %d variants" % (n_dn, i), file=sys.stderr)
            sys.stdout.flush()

        d = denovo(v, samples_lookup, kids, max_alts_in_parents=1, exclude=exclude)
        if d is not None:
            write_denovo(d, sys.stdout)
            n_dn += 1

    print(" called %d de-novos out of %d variants" % (n_dn, i), file=sys.stderr)
    sys.stdout.flush()


def main(argv):
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument("--chrom", help="extract only this chromosome", default=None)
    p.add_argument("--exclude", help="regions to exclude (e.g. LCRs)")
    p.add_argument("--min-mean-depth", default=10, type=int)
    p.add_argument("ped")
    p.add_argument("vcf")

    run(p.parse_args(argv))

if __name__ == "__main__":
    main(sys.argv[1:])
    if __package__ is None:
        from os import path
        sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
