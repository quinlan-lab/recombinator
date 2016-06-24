import argparse
import sys
import re
import itertools as it
import gzip
from collections import defaultdict
from operator import itemgetter, attrgetter

import numpy as np
from cyvcf2 import VCF

from peddy import Ped

HOM_REF, HET, HOM_ALT, UNKNOWN = range(4)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--min-depth", dest='min_depth', type=int, default=20)
    p.add_argument("--min-gq", dest='min_gq', type=int, default=30)
    p.add_argument("--families", default=None, type=str)
    p.add_argument("--ped", required=True)
    p.add_argument("--vcf", required=True)
    p.add_argument("--region", help="optional VCF region e.g. '1:1-1000000'")
    p.add_argument("--prefix", required=True, help="prefix for output. files will be prefix.{region}.{family}.{parent}.bed.gz")
    args = p.parse_args()
    run(args)


def get_family_dict(fam, smp2idx, args):
    """
    Hack to just get the VCF idxs for the dad, mom and kids
    """
    f = {}
    kid_seen = False
    for sample in fam.samples:
        # NOTE: we currently just pull a single quartet from the family...
        if not sample.sex in ("male", "female"): continue
        kids = [s for s in sample.kids if not None in (s.mom, s.dad)]
        if len(kids) < 2: continue

        key = 'dad' if sample.sex == "male" else 'mom'
        f[key] = {'idx': smp2idx[sample.sample_id], 'id': sample.sample_id}

        for i, kid in enumerate(sorted(sample.kids, key=attrgetter('affected'), reverse=True)):
            if i == 2: break # only use first 2 kids
            key = 'sib' if i > 0 else 'template'
            f[key] = {'idx': smp2idx[kid.sample_id], 'id': kid.sample_id}

        other = kid.mom if sample.sex == "male" else kid.dad
        assert other.sample_id != sample.sample_id
        key = 'dad' if other.sex == "male" else 'mom'
        f[key] = {'idx': smp2idx[other.sample_id], 'id': other.sample_id}
        break
    if not f:
        return False

    region = (args.region or "").replace(":", "-")
    f['fh-dad'] = gzip.open("%s.%s.%s.dad.bed.gz" % (args.prefix, region, sample.family_id), "w")
    f['fh-mom'] = gzip.open("%s.%s.%s.mom.bed.gz" % (args.prefix, region, sample.family_id), "w")
    f['fh-dad'].write('\t'.join(['chrom', 'start', 'end', 'parent', 'family_id', 'same', 'dad', 'mom', 'sib1', 'sib2', 'global_call_rate', 'global_depth_1_10_50_90']) + '\n')
    f['fh-mom'].write('\t'.join(['chrom', 'start', 'end', 'parent', 'family_id', 'same', 'dad', 'mom', 'sib1', 'sib2', 'global_call_rate', 'global_depth_1_10_50_90']) + '\n')
    # much faster to index with an array.
    f['idxs'] = np.array([f[s]['idx'] for s in ('dad', 'mom', 'template', 'sib')])
    f['family_id'] = sample.family_id
    return f

def add_genotype_info(fam, gt_types=None,
        gt_depths=None, gt_quals=None):
    """
    Assign the genotype info to each member in the family
    """
    if not gt_types is None:
        fam['gt_type'] = np.array(gt_types[fam['idxs']])
    if not gt_depths is None:
        fam['gt_depth'] = np.array(gt_depths[fam['idxs']])
    if not gt_quals is None:
        fam['gt_qual'] = np.array(gt_quals[fam['idxs']])

def impose_quality_control(fam, args):
    """
    Make sure the genotype data for the family is up to snuff
    """
    # NOTE: any is much faster than np.any for small arrays.
    if any(fam['gt_type'] == UNKNOWN):
        return False

    if any(fam['gt_qual'] < args.min_gq):
        return False

    if any(fam['gt_depth'] < args.min_depth):
        return False

    return True


def is_informative(fam):
    """
    Is the site informative in the sense that to be useful for catching recombination,
    one parent must be HET and the other HOM.
    """
    gt_types = fam['gt_type']  # order is dad, then mom
    return (gt_types[1] == HOM_REF and gt_types[0] == HET) or \
           (gt_types[0] == HOM_REF and gt_types[1] == HET)

def run(args):
    ped = Ped(args.ped)
    vcf = VCF(args.vcf, gts012=True)

    ped_samples = [s.sample_id for s in ped.samples()]
    vcf_samples = set(vcf.samples)

    samples = [s for s in ped_samples if s in vcf_samples]

    vcf = VCF(args.vcf, samples=samples, gts012=True)
    if args.region:
        vcf_iter = vcf(args.region)
    else:
        vcf_iter = vcf

    # build a dict of sample_id to sample index
    smp2idx = dict(zip(vcf.samples, range(len(vcf.samples))))

    # get the Ped objects for the family of interest
    if args.families is None:
        fams = ped.families.values()
    else:
        fams = [ped.families[f] for f in args.families.split(",")]
    if len(fams) == 0:
        sys.exit('Families %s not found in ped file' % args.families)

    # create a simple dictionary of info for each family member
    fs = [get_family_dict(fam, smp2idx, args) for fam in fams]

    # header
    for i, v in enumerate(vcf_iter, start=1):
        if i % 50000 == 0:
            print >>sys.stderr, "at record %d (%s:%d)" % (i, v.CHROM, v.POS)
        if v.var_type != 'snp':
            continue
        if v.call_rate < 0.95: continue

        # expensive to get gt_bases and we only need it at the crossover.
        gt_bases = None
        gt_types, gt_quals, gt_depths = v.gt_types, v.gt_quals, v.gt_depths

        for f in fs:
            # embellish f with the genotype info for each family member.
            # is_informative only needs gt_types, so we check that first...
            add_genotype_info(f, gt_types=gt_types)
            # sanity and quality checks
            if not is_informative(f):
                continue

            # now wee need to add quality and depth.
            add_genotype_info(f, gt_quals=gt_quals, gt_depths=gt_depths)

            if not impose_quality_control(f, args):
                continue

            # detect crossovers.
            for parent, (p1, p2) in [("dad", (0, 1)), ("mom", (1, 0))]:

                if f['gt_type'][p1] == HET and f['gt_type'][p2] == HOM_REF:
                    if gt_bases is None:
                        gt_bases = v.gt_bases
                    fam_bases = "\t".join(gt_bases[f['idxs']])
                    pctiles = "|".join("%.0f" % v for v in
                            np.percentile(v.gt_depths, (1, 10, 50, 90)))

                    val = 1 if f['gt_type'][2] == f['gt_type'][3] else 0
                    f['fh-%s' % parent].write('\t'.join(str(s) for s in [v.CHROM, v.POS - 1, v.POS,
                            parent, f['family_id'], val, fam_bases, "%.2f" %
                            v.call_rate, pctiles]) + '\n')

if __name__ == "__main__":
    main()
