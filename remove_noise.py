import argparse
import sys
import itertools as it
from collections import defaultdict
import toolshed as ts

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--min-size", dest='min_size', type=int, default=50000)
    p.add_argument("--bg", dest='bedgraph', required=True)
    args = p.parse_args()
    run(args)


def report(block, id, keys=['chrom', 'start', 'end', 'state']):
    family_id, parent = id
    print '\t'.join([str(block[k]) for k in keys] + [parent, family_id])


def run(args):
    """
chrom	start	end	parent	family_id	same(1)_diff(2)	dad	mom	sib1	sib2	global_call_rate	global_depth_1_10_50_90
1	10477	10478	mom	11348	2	C/C	C/G	C/G	C/C	1.00	11|24|40|60
1	10488	10489	dad	12290	1	C/T	C/C	C/C	C/C	1.00	10|23|39|55
1	10491	10492	dad	11006	1	C/T	C/C	C/C	C/C	1.00	7|15|33|51
1	10491	10492	mom	11258	1	C/C	C/T	C/T	C/T	1.00	7|15|33|51

    accumulate by state (same) into last. on a state switch, start filling
    next_block. on the next switch, we see if next_block > min_size.
    if yes:
       print last_block and set last_block = next_block
    if no:
       empty next_block and start re-adding to last-block.

    We can do this with using 1 and 2 state only.
    """

    print "\t".join("chrom start end state parent family_id".split())
    last_block = defaultdict(dict)
    next_block = defaultdict(dict)
    block = last_block
    last_chrom = None

    for curr in ts.reader(args.bedgraph):
        if curr['chrom'] == 'chrom' and curr['start'] == 'start': continue
        assert last_chrom is None or last_chrom == curr['chrom'], (last_chrom, curr['chrom'])
        last_chrom = curr['chrom']
        curr['state'] = curr['same(1)_diff(2)']

        id = curr['family_id'], curr['parent']
        if block.get(id) is None:
            block[id]['chrom'] = curr['chrom']
            block[id]['start'] = int(curr['start'])
            block[id]['end'] = int(curr['end'])
            block[id]['state'] = curr['state']
            continue


        prev = block[id]
        if curr['state'] != prev['state']:
            if block == last_block:
                # next_block is empty, so we start filling it.
                block = next_block
                assert block.get(id) is None
                block[id]['chrom'] = curr['chrom']
                block[id]['start'] = int(curr['start'])
                block[id]['end'] = int(curr['end'])
                block[id]['state'] = curr['state']

            else:  # switched states.
                # we check if next_block was big enough so we know if we can
                # drop last_block
                if next_block[id]['end'] - next_block[id]['start'] > args.min_size:
                    # report last_block if it was big enough
                    if last_block[id]['end'] - last_block[id]['start'] > args.min_size:
                        print >>sys.stderr, "reporting last"
                        report(last_block[id], id)

                    last_block = next_block
                # either way, empty next_block
                next_block = defaultdict(dict)
                block = last_block

        else:
            block[id]['end'] = int(curr['end'])

    # cleanup of last blocks
    for block in (next_block, last_block):
        for id, d in block.items():
            if d['end'] - d['start'] > args.min_size:
                report(d, id)

if __name__ == "__main__":
    main()