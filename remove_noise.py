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


def report(block, family_id, parent, keys=['chrom', 'start', 'end', 'state']):
    print '\t'.join([str(block[k]) for k in keys] + [parent, family_id])


def run(args):
    """
chrom	start	end	parent	family_id	same(1)_diff(2)	dad	mom	sib1	sib2	global_call_rate	global_depth_1_10_50_90
1	10477	10478	mom	11348	2	C/C	C/G	C/G	C/C	1.00	11|24|40|60
1	10488	10489	dad	12290	1	C/T	C/C	C/C	C/C	1.00	10|23|39|55
1	10491	10492	dad	11006	1	C/T	C/C	C/C	C/C	1.00	7|15|33|51
1	10491	10492	mom	11258	1	C/C	C/T	C/T	C/T	1.00	7|15|33|51
    """
    last_block = defaultdict(dict)
    all_prev = {}
    curr = {}

    print "\t".join("chrom start end state parent family_id".split())
    for parent in ("mom", "dad"):
        for l in ts.reader(args.bedgraph):
            if l['parent'] != parent: continue
            curr = l
            id = curr['family_id']
            if all_prev.get(id) is None:
                last_block[id]['chrom'] = curr['chrom']
                last_block[id]['start'] = int(curr['start'])
                last_block[id]['end'] = int(curr['end'])
                last_block[id]['state'] = curr['same(1)_diff(2)']
                all_prev[id] = curr
            else:
                prev = all_prev[id]
                if curr['chrom'] != prev['chrom']:
                    report(last_block[id], id, parent)
                else:
                    if curr['same(1)_diff(2)'] != prev['same(1)_diff(2)']:
                        if last_block[id]['end'] - last_block[id]['start'] > args.min_size:
                            report(last_block[id], id, parent)
                        last_block[id] = {}
                        last_block[id]['start'] = int(curr['start'])
                        last_block[id]['chrom'] = curr['chrom']
                        last_block[id]['state'] = curr['same(1)_diff(2)']
                        last_block[id]['end'] = int(curr['end'])
                    else:
                        last_block[id]['end'] = int(curr['end'])
            all_prev[id] = curr
        # cleanup of last blocks
        for id, d in last_block.items():
            report(d, id, parent)

if __name__ == "__main__":
    main()
