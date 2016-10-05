from __future__ import absolute_import, print_function
import sys
import toolshed as ts
from cyvcf2 import VCF
from .var_info import get_position

def main(argv):
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--fragment-length", type=int, default=350, help="library fragment length")
    p.add_argument("--ped-file", required=True, help="optional ped file if parents aren't specified in denovo file")
    p.add_argument("bed", help="bed file of denovos with a column header for 'chrom,start,end,sample_id")
    p.add_argument("vcf")
    args = p.parse_args(argv)
    run(args)

def run(args):

    vcf = VCF(args.vcf, gts012=True)
    sample_lookup = {s: i for i, s in enumerate(vcf.samples)}

    sample_to_dad = {toks[1]: sample_lookup.get(toks[2]) for toks in ts.reader(args.ped_file, header=False)}
    sample_to_mom = {toks[1]: sample_lookup.get(toks[3]) for toks in ts.reader(args.ped_file, header=False)}

    for i, d in enumerate(ts.reader(args.bed, header="ordered")):

        d['mom-sites'] = []
        d['dad-sites'] = []

        idad = sample_to_dad[d['sample_id']]
        imom = sample_to_mom[d['sample_id']]
        ikid = sample_lookup[d['sample_id']]
        for v in get_position(vcf, d, extra=args.fragment_length):
            gt_types = v.gt_types
            if gt_types[ikid] != vcf.HET:
                continue
            if gt_types[idad] == vcf.HET and gt_types[imom] == vcf.HOM_REF:
                d['dad-sites'].append("%s:%d-%d" % (v.CHROM, v.start+1, v.end))
            elif gt_types[imom] == vcf.HET and gt_types[idad] == vcf.HOM_REF:
                d['mom-sites'].append("%s:%d-%d" % (v.CHROM, v.start+1, v.end))

        d['mom-sites'] = ",".join(d['mom-sites'])
        d['dad-sites'] = ",".join(d['dad-sites'])

        if i == 0:
            print("\t".join(d.keys()))
        print("\t".join(d.values()))

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    if __package__ is None:
        from os import path
        sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
    main(sys.argv[1:])
