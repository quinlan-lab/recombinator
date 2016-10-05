from __future__ import absolute_import, print_function
import sys
import toolshed as ts
from cyvcf2 import VCF

def main(argv):
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--fragment-length", type=int, default=300, help="library fragment length")
    p.add_argument("--ped-file", help="optional ped file if parents aren't specified in denovo file")
    p.add_argument("bed", help="bed file of denovos with a column header for 'chrom,start,end,sample_id")
    p.add_argument("vcf")
    args = p.parse_args(argv)
    run(args)

def run(args):

    vcf = VCF(args.vcf, gts012=True)

    if args.ped_file is not None:
        sample_to_dad = {toks[1]: toks[2] for toks in ts.reader(args.ped_file, header=False)}
        sample_to_mom = {toks[1]: toks[3] for toks in ts.reader(args.ped_file, header=False)}
    else:
        sample_to_mom = None
        sample_to_dad = None

    sample_lookup = {s: i for i, s in enumerate(vcf.samples)}

    for i, d in ts.reader(args.bed, header="ordered"):

        d['mom-sites'] = []
        d['dad-sites'] = []

        try:
            idad = sample_lookup[d.get('paternal_id', d.get('dad'))]
        except KeyError:
            idad = sample_to_dad[d['sample_id']]
        try:
            imom = sample_lookup[d.get('maternal_id', d.get('mom'))]
        except KeyError:
            imom = sample_to_mom[d['sample_id']]

        loc = "%s:%d-%d" % (d['chrom'], int(d['start']) + 1 - args.fragment_length, int(d['end'])) + args.fragment_length
        for v in vcf(loc):
            gt_types = v.gt_types
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
