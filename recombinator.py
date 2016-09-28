from __future__ import print_function
import argparse
import time
import shutil
import sys
import os
import re
import itertools as it
import gzip
from collections import defaultdict, OrderedDict
from operator import itemgetter, attrgetter

try:
    import matplotlib.pyplot as plt
    from matplotlib import ticker
    import seaborn as sns
    sns.set_style('whitegrid')
except ImportError:
    pass

import numpy as np
from cyvcf2 import VCF
from intervaltree import IntervalTree

from peddy import Ped

HOM_REF, HET, HOM_ALT, UNKNOWN = range(4)

def filter_main(argv):
    p = argparse.ArgumentParser()
    p.add_argument("--min-sites", type=int, default=20)
    # e.g. we have 2 informatives sites at the supposed break-point, but
    # filtered out a lot of questionable informative sites between them.
    #p.add_argument("--max-intervening", type=int, default=3, help="number of excluded sites inside the actual crossover")
    p.add_argument("--prefix", required=True, help="prefix for output")
    p.add_argument("-p", "--processes", type=int, default=24, help="number of processes")
    p.add_argument("sites", nargs="+", help=".bed.gz files containing state at each informative site")

    args = p.parse_args(argv)
    try:
        os.makedirs(os.path.dirname(args.prefix))
    except OSError:
        pass
    call_all(args.sites, args.prefix, min_sites=args.min_sites,
             processes=args.processes)

def read_exclude(path):
    if path is None:
        return None
    tree = defaultdict(IntervalTree)
    for i, line in enumerate((gzip.open if path.endswith(".gz") else open)(path)):
        toks = line.rstrip().split("\t")
        # skip header if necessary.
        if i == 0:
            try:
                int(toks[1])
            except ValueError:
                continue
        tree[toks[0]].addi(int(toks[1]), int(toks[2]))
    return tree

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--min-depth", dest='min_depth', type=int, default=18)
    p.add_argument("--min-gq", dest='min_gq', type=int, default=20)
    p.add_argument("--families", default=None, type=str)
    p.add_argument("--ped", required=True)
    p.add_argument("--vcf", required=True)
    p.add_argument("--exclude", help="bed file of regions to exclude (e.g. low-complexity)")
    p.add_argument("--region", help="optional VCF region e.g. '1:1-1000000'")
    p.add_argument("--prefix", required=True, help="prefix for output. files will be prefix.{region}.{family}.{parent}.bed.gz")
    args = p.parse_args()

    if args.region:
        args.prefix = os.path.join(args.prefix, args.region.split(":")[0])
    run(args)


def crossovers(f, fhcalls, fhunfilt, prefix, min_sites=20, lock=None):
    """
    Call the actual crossovers. Skip blocks with fewer than min-sites sites and
    then merge the resulting adjacent blocks in the same state.
    """
    fin = gzip.open(f)
    header = next(fin).rstrip().lstrip("#").split("\t")
    ids = [h for h in header if h[-1] == "*"]
    try:
        kid_id = ids[-1].rstrip("*")
        parent_id = ids[0].rstrip("*")
    except IndexError:
        kid_id, parent_id = None, None
    cache = []
    last = None
    for d in (dict(zip(header, l.rstrip().split("\t"))) for l in fin):
        if d['same'] != last:
            if parent_id is None:
                parent_id = d.get('parent_id', None)

            if len(cache) > 0:
                # when we add the next set, we 'collapse' the last block.
                cache[-1] = collapse(cache[-1], parent_id, kid_id)

            cache.append([])

        cache[-1].append(d)
        last = d['same']
    fin.close()

    cache[-1] = collapse(cache[-1], parent_id, kid_id)
    # remove single and double-tons immediately
    cache = remove_bad_regions(cache)
    cache = enforce_min_sites(cache, 1)

    write_crossovers(cache, fhunfilt, lock)
    cache = enforce_min_sites(cache, min_sites)
    xos = write_crossovers(cache, fhcalls, lock)
    xplot(xos, f, prefix)



def find_consecutive_low(vals, min_value=5):
    """
    >>> find_consecutive_low([10, 3, 5])
    []
    >>> find_consecutive_low([10, 2, 2])
    [1, 2]
    >>> find_consecutive_low([10, 2, 2, 8, 2, 2, 8])
    [1, 2, 4, 5]
    >>> find_consecutive_low([10, 2, 2, 8, 2, 2, 2])
    [1, 2, 4, 5, 5, 6]
    >>> find_consecutive_low([2, 2, 2, 8, 2, 2, 2])
    [0, 1, 1, 2, 4, 5, 5, 6]
    >>> find_consecutive_low([2, 2])
    [0, 1]
    >>> find_consecutive_low([2, 10, 2])
    []
    """
    low = []
    for i, (v1, v2) in enumerate(it.izip(vals, vals[1:])):
        if v1 < min_value and v2 < min_value:
            low.append(i)
            low.append(i + 1)
    return low


def remove_bad_regions(cache, n=5):
    """
    remove regions where there is a lot of flipping back and forth
    between states in a short distance.
    """
    assert n % 2 == 1
    for i in range(2):
        if len(cache) < 2: break
        # c0 and c1 contain only values from the same state.
        # if there are adjacent, small blocks, we are pretty sure
        # they are spurious
        c0 = cache[::2]
        assert len(set([c['same'] for c in c0])) == 1, cache
        c1 = cache[1::2]
        assert len(set([c['same'] for c in c1])) == 1, cache

        vals0 = [d['informative-sites'] for d in c0]
        vals1 = [d['informative-sites'] for d in c1]

        drop0 = frozenset(find_consecutive_low(vals0, n))
        drop1 = frozenset(find_consecutive_low(vals1, n))

        cache = [c for i, c in enumerate(c0) if i not in drop0]
        cache += [c for i, c in enumerate(c1) if i not in drop1]

        cache.sort(key=itemgetter('start'))
        # this merges adjacent blocks
        cache = enforce_min_sites(cache, 0)

    #vals = [d['informative-sites'] for d in cache]
    #drop = frozenset(find_consecutive_low(vals, n))
    #cache = [c for i, c in enumerate(cache) if i not in drop]

    return enforce_min_sites(cache, 0)


def write_crossovers(cache, fh, lock=None):
    # now we have starts and ends of the blocks of same state. crossovers occur
    # betwen those blocks.
    if len(cache) < 2: return
    if lock is not None:
        fh = open(fh, mode="a")
        lock.acquire()
    xos = []
    for i in range(len(cache) - 1):
        s, e = cache[i].copy(), cache[i + 1]
        s['left-block'] = "{chrom}:{start}-{end}".format(**s)
        st = s.pop('same')
        s['start'] = int(s['end']) - 1
        s['end'] = e['start']
        s['left-informative-sites'] = s.pop('informative-sites')
        s['right-block'] = "{chrom}:{start}-{end}".format(**e)
        s['right-informative-sites'] = e['informative-sites']
        s['change'] = "%s-%s" % (st, e['same'])
        if i == 0 and fh.tell() == 0:
            # only write the header for the first round.
            print("\t".join(s.keys()), file=fh)
        print("\t".join(map(str, s.values())), file=fh)
        xos.append(s)
    if lock is not None:
        lock.release()
        fh.close()
    return xos


def collapse(dlist, parent_id, kid_id):
    "collapse converts the list of sites into a summary dict of the block."
    assert len(set(d['same'] for d in dlist)) == 1
    d = OrderedDict([
        ('chrom', dlist[0]['chrom']),
        ('start', int(dlist[0]['start'])),
        ('end', int(dlist[-1]['end'])),
        ('same', dlist[0]['same']),
        ('family_id', dlist[0]['family_id']),
        ])

    if parent_id is not None:
        d['parent_id'] = parent_id
    if kid_id is not None:
        d['kid_id'] = kid_id
    d['informative-sites'] = len(dlist)
    return d

def enforce_min_sites(cache, min_sites):
    """Remove blocks with fewer than min_sites informative sites and then
    merge adjacent blocks in the same state."""

    cache = [c for c in cache if c['informative-sites'] >= min_sites]
    if len(cache) < 2: return cache
    icache = [cache[0]]
    for i, c in enumerate(cache[1:]):
        if c['same'] != icache[-1]['same'] or c['chrom'] != icache[-1]['chrom']:
            icache.append(c)
            continue

        icache[-1]['end'] = c['end']
        icache[-1]['informative-sites'] += c['informative-sites']

    return icache

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

    region = (args.region or "all").replace(":", "-")
    try:
        os.makedirs("%s/fam%s" % (args.prefix, sample.family_id))
    except OSError:
        if not os.path.exists("%s/fam%s" % (args.prefix, sample.family_id)):
            raise
    # unphased file-handles
    header = ['#chrom', 'start', 'end', 'parent_id', 'family_id', 'same',
              'dad', 'mom', 'sib1', 'sib2', 'family_depth', 'global_call_rate',
              'global_depth_1_10_50_90', 'family_allele_balance']
    f['fh-dad'] = gzip.open("%s/fam%s/%s.dad.bed.gz" % (args.prefix, sample.family_id, region), "w")
    f['fh-mom'] = gzip.open("%s/fam%s/%s.mom.bed.gz" % (args.prefix, sample.family_id, region), "w")
    f['fh-dad'].write('\t'.join(header) + '\n')
    f['fh-mom'].write('\t'.join(header) + '\n')

    f['ids'] = [f[s]['id'] for s in ('dad', 'mom', 'template', 'sib')]

    # we lable the columns with the 4 sample-ids with a '*' to indicate the
    # parent-child pair for the current scenario.
    for pi, p in enumerate(f['ids'][:2], start=1):
        md = 'dad' if pi == 1 else 'mom'
        dad_lbl = f['ids'][0] + ("*" if pi == 1 else "")
        mom_lbl = f['ids'][1] + ("*" if pi == 2 else "")
        for i, kid in enumerate(f['ids'][2:], start=1):
            sib1_lbl = f['ids'][2] + ("*" if i == 1 else "")
            sib2_lbl = f['ids'][3] + ("*" if i == 2 else "")

            # phased file-handles
            f['fh-%s-%s' % (md, kid)] = gzip.open("%s/fam%s/%s.%s.%s-%s.bed.gz" % (args.prefix, sample.family_id, region, md, p, kid), "w")

            f['fh-%s-%s' % (md, kid)].write('\t'.join(['#chrom', 'start', 'end',
                'parent_id', 'family_id', 'same', dad_lbl, mom_lbl, sib1_lbl, sib2_lbl]) + '\n')

    # much faster to index with an array.
    f['idxs'] = np.array([f[s]['idx'] for s in ('dad', 'mom', 'template', 'sib')])
    f['family_id'] = sample.family_id
    return f

def add_genotype_info(fam, gt_types=None,
        gt_depths=None, gt_quals=None, gt_phases=None):
    """
    Assign the genotype info to each member in the family
    """
    if not gt_phases is None:
        fam['gt_phase'] = gt_phases
    if not gt_types is None:
        fam['gt_type'] = list(gt_types[fam['idxs']])
    if not gt_depths is None:
        fam['gt_depth'] = np.array(gt_depths[fam['idxs']])
    if not gt_quals is None:
        fam['gt_qual'] = np.array(gt_quals[fam['idxs']])

def impose_quality_control(fam, args):
    """
    Make sure the genotype data for the family is up to snuff
    """
    # NOTE: any is much faster than np.any for small arrays.
    if UNKNOWN in fam['gt_type']:
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
    return (gt_types[1] in (HOM_ALT, HOM_REF) and gt_types[0] == HET) or \
           (gt_types[0] in (HOM_REF, HOM_ALT) and gt_types[1] == HET)


def get_allele_balance(v, has_ab):
    if has_ab:
        try:
            sample_abs = v.format("AB")
            return sample_abs
        except KeyError:
            # if "AB" isn't present, use depths directly.
            return None

    ad = (v.gt_alt_depths).astype(float)
    ad[ad < 0] = np.nan
    sample_abs = ad / (ad + v.gt_ref_depths)
    sample_abs[sample_abs < 0] = np.nan
    return sample_abs

def run(args):
    ped = Ped(args.ped)
    vcf = VCF(args.vcf, gts012=True)

    ped_samples = [s.sample_id for s in ped.samples()]
    vcf_samples = set(vcf.samples)

    samples = [s for s in ped_samples if s in vcf_samples]

    exclude = read_exclude(args.exclude)

    vcf = VCF(args.vcf, samples=samples, gts012=True)
    if args.region:
        vcf_iter = vcf(args.region)
    else:
        vcf_iter = vcf

    pctile1 = 10
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
    del fam

    fsites = open("%s.sites" % args.prefix, "w")
    # fcalls contains the crossovers for all samples.
    try:
        vcf["AB"]
        has_abs = True
    except KeyError:
        has_abs = False


    nused, i, report_at, t0 = 0, 0, 10000, time.time()
    for i, v in enumerate(vcf_iter, start=1):
        if i % report_at == 0:
            persec = i / float(time.time() - t0)
            print("%s:%d (%.1f/sec) %.2f%% informative (%d/%d variants)" % (v.CHROM, v.POS,
                  persec, 100.0 * nused/i, nused, i), file=sys.stderr)
            if i == 20000:
                report_at = 40000
            if i == 40000:
                report_at = 100000
            if i == 100000:
                report_at = 200000
                for f in fs:
                    for k in f:
                        if k.startswith('fh'): f[k].flush()
            sys.stderr.flush()
        if v.var_type != 'snp':
            if len(v.REF) > 3 or len(v.ALT) > 1 or len(v.ALT[0]) > 3:
                continue
        if v.call_rate < 0.95: continue
        if v.FILTER is not None: continue

        if exclude is not None and 0 != len(exclude[v.CHROM].search(v.start, v.end)):
            continue

        # expensive to get gt_bases and we only need it at the crossover.
        gt_bases = None
        gt_types, gt_quals, gt_depths = v.gt_types, v.gt_quals, v.gt_depths
        gt_phases = v.gt_phases
        ipctiles, pctiles = None, None
        sample_abs = None

        nsites = 0 # track the number of families that had this as an informative site.
        for f in fs:
            if ipctiles is not None and ipctiles[0] < pctile1:
                break

            # is_informative only needs gt_types, so we check that first...
            add_genotype_info(f, gt_types=gt_types, gt_phases=gt_phases)

            # ############## PHASED ####################
            if all(f['gt_phase']):
                if gt_bases is None:
                    gt_bases = v.gt_bases

                nsites += phased_check(f, v, gt_bases)
                continue
            # ############ END PHASED ####################

            # need exactly 1 het parent for unphased checks.
            if 1 != ((f['gt_type'][0] == HET) + (f['gt_type'][1] == HET)):
                continue

            if not is_informative(f):
                continue

            # now wee need to add quality and depth.
            add_genotype_info(f, gt_quals=gt_quals, gt_depths=gt_depths)

            if not impose_quality_control(f, args):
                continue

            # detect crossovers.
            for parent, (p1, p2) in [("dad", (0, 1)), ("mom", (1, 0))]:

                if not (f['gt_type'][p1] == HET and f['gt_type'][p2] == HOM_REF):
                    continue

                if gt_bases is None:
                    gt_bases = v.gt_bases
                if sample_abs is None:
                    sample_abs = get_allele_balance(v, has_abs)
                    if sample_abs is None: break

                fam_abs = sample_abs[f['idxs']]
                off = 0.31  # require that  off <= alt/(ref+alt) <= 1-off
                if ((fam_abs[p1] >= 1 - off) | (fam_abs[p1] <= off)): continue
                if np.any((1 - off < fam_abs[2:]) | (fam_abs[2:] <= off)): continue


                fam_bases = "\t".join(gt_bases[f['idxs']])

                fam_abs = "|".join("%.2f" % val for val in fam_abs)

                # calculate on first use. we found that having a low 1st pctile
                # was a good indicator of increased chance of spurious XO even
                # in families with decent depth.
                if pctiles is None:
                    ipctiles = np.percentile(gt_depths, (1, 10, 50, 90))
                    pctiles = "|".join("%.0f" % de for de in ipctiles)
                if ipctiles[0] < pctile1:
                    break

                fam_depths = "|".join(map(str, gt_depths[f['idxs']]))
                nsites += 1
                val = 1 if f['gt_type'][2] == f['gt_type'][3] else 0
                f['fh-%s' % parent].write('\t'.join(str(s) for s in [v.CHROM, v.POS - 1, v.POS,
                        f['ids'][p1], f['family_id'], val, fam_bases,
                        fam_depths, "%.2f" % v.call_rate, pctiles, fam_abs]) + '\n')

        fsites.write("%s:%d\t%d\n" % (v.CHROM, v.POS, nsites))
        if nsites > 0:
            nused += 1

    fsites.close()
    persec = i / float(time.time() - t0)
    print("finished at %s:%d (%.1f/sec) %.2f%% informative (%d/%d variants)" % (v.CHROM, v.POS,
        persec, 100.0 * nused/i, nused, i), file=sys.stderr)
    kept = _remove_empty(fs)
    call_all(kept, args.prefix, min_sites=20)


def call_all(kept, prefix, min_sites=20, processes=1):
    """
    call all takes the informative sites, calls the crossovers, and makes plots.
    """
    iprefix = prefix if prefix[-1] in "./" else (prefix + ".")

    fcalls = open("%scrossovers.bed" % iprefix, "w")
    funfiltered = open("%scrossovers-unfiltered.bed" % iprefix, "w")

    lock, pool = None, None
    if processes > 1:
        from concurrent import futures
        # note: may want to switch to a different file-lock, but this seems
        # fine.
        # http://fasteners.readthedocs.io/en/latest/api/process_lock.html#classes
        import lockfile
        pool = futures.ProcessPoolExecutor(processes)
        fcalls.close()
        funfiltered.close()
        fcalls = fcalls.name
        funfiltered = funfiltered.name
        lock = lockfile.LockFile(fcalls + ".lock")

    t0, n = time.time(), 1000
    for i, f in enumerate(kept, start=1):
        if i % n == 0:
            persec = n / float(time.time() - t0)
            t0 = time.time()
            print("%d/%d (%.3f/sec)" % (i, len(kept), persec), file=sys.stderr)
        if processes == 1:
            crossovers(f, fcalls, funfiltered, prefix, min_sites, lock)
        else:
            pool.submit(crossovers, f, fcalls, funfiltered, prefix, min_sites, lock)

    if pool is not None:
        pool.shutdown()
        print("wrote aggregated calls to: %s and %s" % (fcalls, funfiltered), file=sys.stderr)
    else:
        fcalls.close()
        funfiltered.close()
        print("wrote aggregated calls to: %s and %s" % (fcalls.name, funfiltered.name), file=sys.stderr)


def rdr(f, ordered=False):

    fh = (gzip.open if f.endswith(".gz") else open)(f)
    header = next(fh).rstrip("\r\n").split("\t")
    if header[0][0] == "#":
        header[0] = header[0][1:]

    d = OrderedDict if ordered else dict

    for toks in (l.rstrip("\r\n").split("\t") for l in fh):
        yield d(zip(header, toks))
    fh.close()


def xplot(xos, fsites, prefix):
    figname = os.path.basename(fsites).replace(".bed.gz", ".png")
    row = next(rdr(fsites))
    if prefix.endswith(os.path.sep + row['chrom']):
        d = os.path.sep.join((prefix, "fam" + row['family_id']))
    else:
        d = os.path.sep.join((prefix, row['chrom'], "fam" + row['family_id']))
    try:
        os.makedirs(d)
    except OSError:
        pass
    figname = os.path.join(d, figname)

    fig, ax = plt.subplots(1, sharex=True, figsize=(7, 2))

    name = fsites.split("/")[-1].rsplit(".", 2)[0]

    xs, ys = [], []
    for row in rdr(fsites):
        xs.append(int(row['start']))
        ys.append(int(row['same']) / 2.0 + 0.25) # 0 -> 0.25, 1 -> 0.75

    if name in ('mom', 'dad') or name.endswith((".mom", ".dad")):
        name += " " + row['parent_id']

    if not row['family_id'] in name:
        name += " (family: %s)" % row['family_id']
    ax.set_title(name)
    ax.plot(xs, ys, color='k', ls='none', marker='.', markersize=5, zorder=5)
    rng = xs[-1] - xs[0]

    def fmtr(x, p):
        if rng > 3000000:
            v, suf = "%.2f" % (x / 1000000.), "M"
        else:
            v, suf = "%.2f" % (x / 1000.), "K"
        # strip trailing 0's and "."
        while v[-1] in '.0' and '.' in v:
            v = v[:-1]
        return v + suf

    ax.get_xaxis().set_major_formatter(ticker.FuncFormatter(fmtr))
    ax.set_xlim(xmin=xs[0] - 1, xmax=xs[-1] + 1)
    ax.set_ylim(ymin=0, ymax=1)
    ax.set_xlabel('Genomic Position')
    ax.set_yticks([0.25, 0.75])
    ax.set_yticklabels(["state 0", "state 1"], rotation='vertical')

    if xos is None:
        plt.tight_layout()
        plt.savefig(figname)
        plt.close(fig)
        return

    left = xs[0]
    # plot the runs and plot the actual crossovers in red.
    for row in xos:
        left_state = int(row['change'].split("-")[0])
        ax.axvspan(xmin=left, xmax=int(row['start']),
                ymin=left_state / 2., ymax=(left_state+1) / 2.,
                ec='none')
        left = int(row['end'])


        # plot the actual crossover:
        ax.axvspan(xmin=int(row['start']), xmax=int(row['end']),
                color='#ff3333', zorder=4)

        # plot some extra so we're sure we can see it.
        ax.axvspan(xmin=int(row['start']) - 0.005 * rng, xmax=int(row['end']) +
                0.005 * rng,
                color='#ff3333', zorder=4, alpha=0.3)

    # plot the last block
    right_state = int(row['change'].split("-")[1])
    ax.axvspan(xmin=left, xmax=xs[-1],
            ymin=right_state / 2., ymax=(right_state+1) / 2., ec='none')

    plt.tight_layout()
    plt.savefig(figname)
    plt.close(fig)


def _remove_empty(fs):
    # remove empty files.
    kept = []
    for f in fs:
        for v in f.values():
            if hasattr(v, "flush"):
                v.flush()
                v.close()
                keep = True
                with gzip.open(v.name) as gfh:
                    for i, line in enumerate(gfh):
                        if i == 2:
                            kept.append(v.name)
                            break
                    else:
                        keep = False
                if not keep:
                    os.unlink(v.name)
    return kept


def phased_check(fam, v, gt_bases):
    """
    take cases where only 1 parent is HET
    """
    gt_type = fam['gt_type']
    # need at least 1 het kid
    if 0 == ((gt_type[2] == 1) + (gt_type[3] == 1)):
        return 0

    used = 0
    for parent, (p1, p2) in [("dad", (0, 1)), ("mom", (1, 0))]:
        if gt_type[p1] != HET: continue
        # if gt_type[p2] == HET: continue

        fam_bases = [x.split("|") for x in gt_bases[fam['idxs']]]
        vbases = "\t".join(gt_bases[fam['idxs']])

        for kid in (2, 3):

            if gt_type[kid] != HET: continue
            ref = v.ALT[0]

            try:
                kidx = fam_bases[kid].index(ref)
                pidx = fam_bases[p1].index(ref)
            except ValueError:  # one of them had '.'
                continue

            same = int(kidx == pidx)
            used += 1

            kid_id = fam['ids'][kid]
            parent_id = fam['ids'][p1]
            fam['fh-%s-%s' % (parent, kid_id)].write('\t'.join(str(s) for s in [v.CHROM, v.POS - 1, v.POS,
                    parent_id, fam['family_id'], same, vbases
                    #, "%.2f" % v.call_rate, pctiles
                    ]) + '\n')
    return int(used > 0)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    if len(sys.argv) > 1 and sys.argv[1] == "filter":
        sys.exit(filter_main(sys.argv[2:]))

    main()
