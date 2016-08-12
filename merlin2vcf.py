import itertools as it

def tryint(f):
    try:
        return str(int(f) - 1)
    except:
        if "?" == f:
            return "?"
        if "," in f:
            return "HET"
        else:
            raise

def _get_fam(fh):
    try:
        fam = next(fh)
    except StopIteration:
        return None, None
    fam = fam[7:].split("[", 1)[0].strip()
    samples = []
    line = next(fh).rstrip()
    while line:
        line = line.split()
        line[1] = line[0]
        line = line[1:]
        line[3:] = map(tryint, line[3:])
        samples.append(line)
        line = next(fh).rstrip()
    return fam, samples


def _read_map(f):
    return [x[1].split("_") for x in (l.rstrip().split() for i, l in enumerate(open(f)) if i > 0)]


def merlin2vcf(fh, fmap):
    L = []
    fam, samps = _get_fam(fh)
    while fam:
        L.extend(samps)
        fam, samps = _get_fam(fh)

    pmap = _read_map(fmap)
    samples = "\t".join(set([s[0] for s in L]))
    print """\
##fileformat=VCFv4.1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	{samples}""".format(samples=samples)
    # start at 1 to skip the sample.
    for i, p in enumerate(pmap, start=1):
        try:
            gts = [L[k][i] + "|" + L[k+1][i] for k in range(0, len(L), 2)]
        except:
            print i, len(s)
            raise
        for k, gt in enumerate(gts):
            if "?" in gt or "," in gt:
                gts[k] = "./."
            elif "HET" in gt:
                gts[k] = "0/1"

        print "%s\t%s\t.\tA\tC\t20\tPASS\t.\tGT\t%s" % (p[0], p[1], "\t".join(gts))


if __name__ == "__main__":
    import sys
    merlin2vcf(open(sys.argv[1]), sys.argv[2])
