from __future__ import print_function
import sys
import toolshed as ts
import itertools as it
from operator import itemgetter
from concurrent.futures import ProcessPoolExecutor

fam_beds = sys.argv[1:]

min_sites = 2
min_flank = 20
max_sites = 20
max_length = 3000


print("#chrom\tstart\tend\tsites-left\tsites-in\tsites-right\tparent_id\tfamily_id\tab\tdepth\tfile")

def get_ab(grp):
    ab = []
    for d in grp:
        ab.append("|".join(x for x in d['family_allele_balance'].split("|") if x != 'nan'))
    return ",".join(ab)

def find_gcs(fam_bed):
    last_grp = None
    ends = []
    nsites = []
    res = []

    for _, grp in it.groupby(ts.reader(fam_bed), itemgetter('same')):

        grp = list(grp)
        ends.append(int(grp[-1]['end']))
        nsites.append(len(grp))
        # cant find until we have at least 3 groups.
        if len(ends) < 3:
            last_grp = grp
            continue

        start = int(grp[0]['start'])
        # check num sites, distance between flanking regions and number of
        if min_sites <= len(last_grp) <= max_sites and (start - ends[-3]) < max_length and nsites[-1] >= min_flank and nsites[-3] >= min_flank:
            d = last_grp[0]

            res.append("\t".join(map(str, (d['chrom'], ends[-3], start, nsites[-3],
                    nsites[-2], nsites[-1], d['parent_id'], d['family_id'],
                    get_ab(last_grp),
                    ",".join(x['family_depth'] for x in last_grp),
                    fam_bed))))
            assert nsites[-2] == len(last_grp)

        # set it here, but still need to check if it's bounds are small enough.
        last_grp = grp

    return res

with ProcessPoolExecutor(10) as p:
    for lines in p.map(find_gcs, (f for f in fam_beds)):
        if len(lines) > 0:
            print("\n".join(lines))
            sys.stdout.flush()
