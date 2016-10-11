import sys
import toolshed as ts
import itertools as it
from operator import itemgetter
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Manager

fam_beds = sys.argv[1:]

min_sites = 2
min_flank = 10
max_sites = 20
max_length = 3000


print("#chrom\tstart\tend\tsites-left\tsites-in\tsites-right\tparent_id\tfamily_id")

def find_gcs(args):
    fam_bed, lock = args
    last_grp = None
    ends = []
    nsites = []
    try:
        for _, grp in it.groupby(ts.reader(fam_bed, header="ordered"), itemgetter('same')):

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
                d = grp[0]

                ab_fail = False
                for d in grp:
                    if any(not 0.40 < float(v) < 0.60 for v in d['family_allele_balance'].split("|") if v != "nan"):
                        ab_fail = True
                        break

                if not ab_fail:
                    lock.acquire()
                    print("\t".join(map(str, (d['chrom'], ends[-3], start, nsites[-3],
                        nsites[-2], nsites[-1], d['parent_id'], d['family_id'], fam_bed))))
                    lock.release()
                assert nsites[-2] == len(last_grp)

            # set it here, but still need to check if it's bounds are small enough.
            last_grp = grp
    except:
        pass

with ProcessPoolExecutor() as p:
    l = Manager().Lock()
    list(p.map(find_gcs, ((f, l) for f in fam_beds)))
