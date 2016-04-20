import argparse
import sys
import re
import itertools as it
from collections import defaultdict
from cyvcf2 import VCF

from pedagree import Ped

HOM_REF = 0
HET = 1
HOM_ALT = 3
UNKNOWN = 2


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--min-depth", dest='min_depth', type=int, default=20)
    p.add_argument("--min-gq", dest='min_gq', type=int, default=30)
    p.add_argument("--family", dest='family', required=True)
    p.add_argument("--parent", dest='parent',
                   default="dad", choices=['dad', 'mom'])
    p.add_argument("--ped", dest='ped', required=True)
    p.add_argument("--vcf", dest='vcf', required=True)
    args = p.parse_args()
    run(args)


def get_family_dict(fam, smp2idx):
    """
    Hack to just get the VCF idxs for the dad, mom and kids
    """
    f = defaultdict(defaultdict)
    kid_seen = False
    for sample in fam.samples:
        if sample.dad is None and sample.mom is None:
            if sample.sex == "male":
                f['dad']['idx'] = smp2idx[sample.sample_id]
                f['dad']['id'] = sample.sample_id
            elif sample.sex == "female":
                f['mom']['idx'] = smp2idx[sample.sample_id]
                f['mom']['id'] = sample.sample_id
        else:
            # first kid is template.
            # better? smarter?
            if kid_seen == False:
                f['template']['idx'] = smp2idx[sample.sample_id]
                f['template']['id'] = sample.sample_id
                kid_seen = True
            else:
                f['sib']['idx'] = smp2idx[sample.sample_id]
                f['sib']['id'] = sample.sample_id
    return f


def add_genotype_info(fam, variant):
    """
    Assign the genotype info to each member in the family
    """
    for f in fam:
        fam[f]['gt_type'] = variant.gt_types[fam[f]['idx']]
        fam[f]['gt_base'] = variant.gt_bases[fam[f]['idx']]
        fam[f]['gt_phase'] = variant.gt_phases[fam[f]['idx']]
        fam[f]['gt_depth'] = variant.gt_depths[fam[f]['idx']]
        fam[f]['gt_qual'] = variant.gt_quals[fam[f]['idx']]


def impose_quality_control(fam, args):
    """
    Make sure the genotype data for the family is up to snuff
    """
    gt_types = [fam[f]['gt_type'] for f in fam]
    gt_quals = [fam[f]['gt_qual'] for f in fam]
    gt_depths = [fam[f]['gt_depth'] for f in fam]

    if any(g == UNKNOWN for g in gt_types):
        return False
    if any(g < args.min_gq for g in gt_quals):
        return False
    if any(g < args.min_depth for g in gt_depths):
        return False
    return True


def is_informative(fam):
    """
    Is the site informative in the sense that
    to be useful for catching recombination,
    one parent must be HET and the other HOM
    """
    if (fam['mom']['gt_type'] == HOM_REF and fam['dad']['gt_type'] == HET) \
            or \
            (fam['mom']['gt_type'] == HET and fam['dad']['gt_type'] == HOM_REF):
        return True
    else:
        return False


def run(args):
    vcf = VCF(args.vcf)
    ped = Ped(args.ped)

    # build a dict of sample_id to sample index
    smp2idx = dict(zip(vcf.samples, range(len(vcf.samples))))

    # get the Ped objects for the family of interest
    fam = ped.families.get(args.family)
    if fam is None:
        sys.exit('Family %s not found in ped file' % args.family)

    # create a simple dictionary of info for each family member
    f = get_family_dict(fam, smp2idx)

    # header
    print '\t'.join(['chrom', 'start', 'end', 'same(1)_diff(2)',
                     'dad_' + f['dad']['id'], 'mom_' + f['mom']['id'],
                     'template_' + f['template']['id'], 'sib_' + f['sib']['id']])
    for v in vcf:
        # embellish f with the genotype info for each family member.
        add_genotype_info(f, v)

        # sanity and quality checks
        if not v.var_type == 'snp':
            continue
        if not impose_quality_control(f, args):
            continue
        if not is_informative(f):
            continue

        # detect crossovers.
        p1 = "dad"
        p2 = "mom"
        if args.parent == 'mom':
            p1, p2 = p2, p1

        if f[p1]['gt_type'] == HET and f[p2]['gt_type'] == HOM_REF:
            if f['template']['gt_type'] == f['sib']['gt_type']:
                print '\t'.join(str(s) for s in [v.CHROM, v.POS - 1, v.POS, 1, f['dad']['gt_base'], f['mom']['gt_base'], f['template']['gt_base'], f['sib']['gt_base']])
            else:
                print '\t'.join(str(s) for s in [v.CHROM, v.POS - 1, v.POS, 2, f['dad']['gt_base'], f['mom']['gt_base'], f['template']['gt_base'], f['sib']['gt_base']])

if __name__ == "__main__":
    main()
