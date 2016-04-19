import argparse
import re
import itertools as it
from collections import defaultdict
from cyvcf2 import VCF
from pedagree import Ped


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--min-depth", dest='min_depth', type=int, default=20)
    p.add_argument("--min-gq", dest='min_gq', type=int, default=30)
    p.add_argument("--family", dest='family', required=True)
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
            elif sample.sex == "female":
                f['mom']['idx'] = smp2idx[sample.sample_id]
        else:
            # first kid is template.
            # better? smarter?
            if kid_seen == False:
                f['template']['idx'] = smp2idx[sample.sample_id]
                kid_seen = True
            else:
                f['sib']['idx'] = smp2idx[sample.sample_id]
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

    if any(g == 2 for g in gt_types):
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
    if (fam['mom']['gt_type'] == 0 and fam['dad']['gt_type'] == 1) \
         or \
         (fam['mom']['gt_type'] == 1 and fam['dad']['gt_type'] == 0):
        return True
    else:
        return False


def run(args):
    vcf = VCF(args.vcf)
    ped = Ped(args.ped)

    # build a dict of sample_id to sample index
    smp2idx = dict(zip(vcf.samples, range(len(vcf.samples))))

    # get the Ped objects for the family of interest
    fam = ped.families[args.family]

    # create a simple dictionary of info for each family member
    f = get_family_dict(fam, smp2idx)

    for v in vcf:
        # embellish f with the genotype info for each family member.
        add_genotype_info(f, v)
    
        # sanity and quality checks
        if not (len(v.REF) == 1 and len(v.ALT) == 1):
            continue
        if not impose_quality_control(f, args):
            continue
        if not is_informative(f):
            continue

        # detect paternal crossovers.
        # TODO: maternal.
        if f['dad']['gt_type'] == 1 and f['mom']['gt_type'] == 0:
            if f['template']['gt_type'] == f['sib']['gt_type']:
                print '\t'.join(str(s) for s in [v.CHROM, v.POS-1, v.POS, 1, f['dad']['gt_base'], f['mom']['gt_base'], f['template']['gt_base'], f['sib']['gt_base']])
            else:
                print '\t'.join(str(s) for s in [v.CHROM, v.POS-1, v.POS, 2, f['dad']['gt_base'], f['mom']['gt_base'], f['template']['gt_base'], f['sib']['gt_base']])

if __name__ == "__main__":
    main()