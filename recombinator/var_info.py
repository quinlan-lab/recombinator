from __future__ import absolute_import, print_function
import sys
import toolshed as ts
from .denovo import variant_info
from peddy import Ped
from cyvcf2 import VCF

def main(argv):
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("bed", help="bed file of variants for which to extract info")
    p.add_argument("ped")
    p.add_argument("vcf")
    args = p.parse_args(argv)
    run(args)

def get_position(vcf, d, ref=None, alt=None, extra=0):
    loc = "%s:%d-%d" % (d['chrom'], int(d['start']) + 1 - extra,
                        int(d['end']) + extra)
    for v in vcf(loc):
        if ref and v.REF != ref: continue
        if alt and alt not in v.ALT: continue
        yield v

def run(args):

    vcf = VCF(args.vcf, gts012=True)
    sample_lookup = {s: i for i, s in enumerate(vcf.samples)}

    ped = Ped(args.ped)

    kid_lookup = {s.sample_id: s for s in ped.samples()}

    for i, d in enumerate(ts.reader(args.bed, header="ordered")):

        rec = None
        kid = kid_lookup[d['sample_id']]
        for v in get_position(vcf, d):
            rec = variant_info(v, kid, sample_lookup)

        if i == 0:
            print("\t".join(rec.keys()))
        print("\t".join(map(str, rec.values())))

if __name__ == "__main__":
    if __package__ is None:
        from os import path
        sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
    import doctest
    doctest.testmod()
    main(sys.argv[1:])
