from __future__ import print_function
import argparse
import sys
import os
import re
import itertools as it
import gzip
from collections import defaultdict, OrderedDict
from operator import itemgetter, attrgetter

import numpy as np
from cyvcf2 import VCF

from peddy import Ped

HOM_REF, HET, HOM_ALT, UNKNOWN = range(4)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--min-depth", dest='min_depth', type=int, default=18)
    p.add_argument("--min-gq", dest='min_gq', type=int, default=20)
    p.add_argument("--families", default=None, type=str)
    p.add_argument("--ped", required=True)
    p.add_argument("--vcf", required=True)
    p.add_argument("--region", help="optional VCF region e.g. '1:1-1000000'")
    p.add_argument("--prefix", required=True, help="prefix for output. files will be prefix.{region}.{family}.{parent}.bed.gz")
    args = p.parse_args()

    if args.region:
        args.prefix = os.path.join(args.prefix, args.region.split(":")[0])
    run(args)

def crossovers(f, min_sites=20):
    """
    call the actual crossovers. skip blocks with fewer than min-sites sites and
    then merge the resulting adjacent blocks in the same state
    """
    fin = gzip.open(f)
    header = next(fin).rstrip().split("\t")
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
                cache[-1] = collapse(cache[-1], parent_id, kid_id)

            cache.append([])

        cache[-1].append(d)
        last = d['same']
    fin.close()

    cache[-1] = collapse(cache[-1], parent_id, kid_id)
    # remove single and double-tons immediately
    cache = remove_bad_regions(cache)
    cache = enforce_min_sites(cache, 3)
    write_crossovers(cache, f.replace('.bed.gz', '.crossovers-unfiltered.bed'))
    cache = enforce_min_sites(cache, min_sites)
    return write_crossovers(cache, f.replace('.bed.gz', '.crossovers.bed'))


def xmean(x, N):
    """
    >>> xmean(np.array([1, 2, 3, 2, 1, 600, 2, 10]), 2)
    array([   2.        ,    2.33333333,    2.        ,  201.        ,
            201.        ,  204.        ,    6.        ,   10.        ])
    """

    res = []
    for i in range(len(x)):
        res.append(np.mean(x[i: i+N+1]))
    res = np.array(res)
    assert res.shape == x.shape
    return res


def remove_bad_regions(cache, n=5):
    """
    remove regions where there is a lot of flipping back and forth
    between states in a short distance.
    """
    assert n % 2 == 1
    for i in range(2):
        if len(cache) < 2: break
        c0 = cache[::2]
        assert len(set([c['same'] for c in c0])) == 1, cache
        c1 = cache[1::2]
        assert len(set([c['same'] for c in c1])) == 1, cache

        vals0 = np.array([d['informative-sites'] for d in c0])
        vals1 = np.array([d['informative-sites'] for d in c1])

        keep0, = np.where(xmean(vals0, 2) >= n)
        keep1, = np.where(xmean(vals1, 2) >= n)
        keep0, keep1 = set(keep0), set(keep1)

        cache = [c for i, c in enumerate(c0) if i in keep0]
        cache += [c for i, c in enumerate(c1) if i in keep1]

        cache.sort(key=itemgetter('start'))
        # this merges adjacent blocks
        cache = enforce_min_sites(cache, 0)
    return enforce_min_sites(cache, 0)

def write_crossovers(cache, fname):
    # now we have starts and ends of the blocks of same state. crossovers occur
    # betwen those blocks.
    if len(cache) < 2: return
    o = open(fname, "w")
    for i in range(len(cache) - 1):
        s, e = cache[i].copy(), cache[i + 1]
        st = s.pop('same')
        s['start'] = int(s['end']) - 1
        s['end'] = e['start']
        s['informative-sites-r'] = e['informative-sites']
        s['change'] = "%s-%s" % (st, e['same'])
        if i == 0:
            print("\t".join(s.keys()), file=o)
        print("\t".join(map(str, s.values())), file=o)
    o.close()
    return o.name


def collapse(dlist, parent_id, kid_id):
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
    f['fh-dad'] = gzip.open("%s/fam%s/%s.dad.bed.gz" % (args.prefix, sample.family_id, region), "w")
    f['fh-mom'] = gzip.open("%s/fam%s/%s.mom.bed.gz" % (args.prefix, sample.family_id, region), "w")
    f['fh-dad'].write('\t'.join(['chrom', 'start', 'end', 'parent_id', 'family_id', 'same', 'dad', 'mom', 'sib1', 'sib2', 'global_call_rate', 'global_depth_1_10_50_90']) + '\n')
    f['fh-mom'].write('\t'.join(['chrom', 'start', 'end', 'parent_id', 'family_id', 'same', 'dad', 'mom', 'sib1', 'sib2', 'global_call_rate', 'global_depth_1_10_50_90']) + '\n')

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

            f['fh-%s-%s' % (md, kid)] = gzip.open("%s/fam%s/%s.%s.%s-%s.bed.gz" % (args.prefix, sample.family_id, region, md, p, kid), "w")

            f['fh-%s-%s' % (md, kid)].write('\t'.join(['chrom', 'start', 'end',
                'family_id', 'same', dad_lbl, mom_lbl, sib1_lbl, sib2_lbl,
                #'global_call_rate', 'global_depth_1_10_50_90'
                ]) + '\n')

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
    return (gt_types[1] in (HOM_ALT, HOM_REF) and gt_types[0] == HET) or \
           (gt_types[0] in (HOM_REF, HOM_ALT) and gt_types[1] == HET)

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
    del fam

    # header
    for i, v in enumerate(vcf_iter, start=1):
        if i % 200000 == 0:
            print("at record %d (%s:%d)" % (i, v.CHROM, v.POS), file=sys.stderr)
        if v.var_type != 'snp':
            if len(v.REF) > 3 or len(v.ALT) > 1 or len(v.ALT[0]) > 3:
                continue
        if v.call_rate < 0.95: continue

        # expensive to get gt_bases and we only need it at the crossover.
        gt_bases = None
        gt_types, gt_quals, gt_depths = v.gt_types, v.gt_quals, v.gt_depths
        gt_phases = v.gt_phases

        for f in fs:
            # embellish f with the genotype info for each family member.
            # is_informative only needs gt_types, so we check that first...
            add_genotype_info(f, gt_types=gt_types, gt_phases=gt_phases)
            # sanity and quality checks

            if np.all(f['gt_phase']):
                if gt_bases is None:
                    gt_bases = v.gt_bases

                phased_check(f, v, gt_bases)
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
                fam_bases = "\t".join(gt_bases[f['idxs']])
                pctiles = "|".join("%.0f" % v for v in
                        np.percentile(v.gt_depths, (1, 10, 50, 90)))

                val = 1 if f['gt_type'][2] == f['gt_type'][3] else 0
                f['fh-%s' % parent].write('\t'.join(str(s) for s in [v.CHROM, v.POS - 1, v.POS,
                        f['ids'][p1], f['family_id'], val, fam_bases, "%.2f" %
                        v.call_rate, pctiles]) + '\n')

    kept = _remove_empty(fs)
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
        plot = xplot
    except ImportError:
        sys.stderr.write("install matplotlib and seaborn for plots\n")
        plot = _plot
    for f in kept:
        xos = crossovers(f)
        plot(xos, f)

def _plot(*args, **kwargs):
    pass

def rdr(f, ordered=False):

    fh = (gzip.open if f.endswith(".gz") else open)(f)
    header = next(fh).rstrip("\r\n").split("\t")
    if header[0][0] == "#":
        header[0] = header[0][1:]

    for toks in (l.rstrip("\r\n").split("\t") for l in fh):
        if ordered:
            yield OrderedDict(zip(header, toks))
        else:
            yield dict(zip(header, toks))
    fh.close()

def xplot(fxos, fsites):
    import matplotlib.pyplot as plt
    from matplotlib import ticker
    import seaborn as sns
    figname = fsites.replace(".bed.gz", ".png")

    fig, ax = plt.subplots(1, sharex=True)

    ax.set_title("crossovers")
    name = fsites.split("/")[-1].split(".")[0]
    fig.suptitle(name)

    xs, ys = [], []
    for row in rdr(fsites):
        xs.append(int(row['start']))
        ys.append(int(row['same']))

    ax.plot(xs, ys, 'k.', zorder=5)
    rng = xs[-1] - xs[0]

    def fmtr(x, p):
        if xs[-1] - xs[0] > 3000000:
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
    ax.set_ylabel('State')
    ax.set_xlabel('Genomic Position')

    if fxos is None:
        plt.savefig(figname)
        return

    left = xs[0]
    # plot the runs and plot the actual crossovers in red.
    for row in rdr(fxos):
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

    plt.savefig(figname)

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

    for parent, (p1, p2) in [("dad", (0, 1)), ("mom", (1, 0))]:
        if fam['gt_type'][p1] != HET: continue
        if fam['gt_type'][p2] == HET: continue
        # TODO: add impose_quality_control here.

        fam_bases = [x.split("|") for x in gt_bases[fam['idxs']]]
        vbases = "\t".join(gt_bases[fam['idxs']])
        for kid in (2, 3):

            if fam['gt_type'][kid] != HET: continue
            ref = v.ALT[0]

            try:
                kidx = fam_bases[kid].index(ref)
                pidx = fam_bases[p1].index(ref)
            except ValueError:  # one of them had '.'
                continue

            same = int(kidx == pidx)

            kid_id = fam['ids'][kid]
            parent_id = fam['ids'][p1]
            fam['fh-%s-%s' % (parent, kid_id)].write('\t'.join(str(s) for s in [v.CHROM, v.POS - 1, v.POS,
                    fam['family_id'], same, vbases
                    #, "%.2f" % v.call_rate, pctiles
                    ]) + '\n')
            fam['fh-%s-%s' % (parent, kid_id)].flush()

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    main()
