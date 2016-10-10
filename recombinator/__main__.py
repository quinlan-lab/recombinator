from __future__ import print_function, absolute_import
import warnings
import sys

warnings.simplefilter('ignore')
from . import recombinator
from . import cohort_plots
from . import denovo
from . import extract_hq
from . import enrichment
from . import var_info
from . import var_plot
from . import rbp
from . import __version__

import sys

commands = [
        ('recombinator', 'find informative sites and crossovers from a VCF', recombinator.main),
        ('filter', 'call crossovers from informatives sites', recombinator.filter_main),
        ('cohort-plots', 'make plots summarizing crossovers across a cohort', cohort_plots.main),
        ('extract-hq', 'pull high-quality sites (and de-novos) from a vcf to be phased', extract_hq.main),
        ('denovos', 'call de novo variants in the cohort', denovo.main),
        ('enrichment', 'test for enrichment by shuffling labels', enrichment.main),
        ('rbp', 'pull out sites that are candidates for read-backed phasing', rbp.main),
        ('variant-info', 'extract information for all variants in a bed file', var_info.main),
        ('variant-plot', 'plot distribution of values for variation-fino', var_plot.main),
        ]

def print_commands():
    print("recombinator version: %s" % __version__)
    print("""
commands:
%s""" % "\n".join(("%-14s: %s" % (c[0], c[1])) for c in commands))


def main():
    if len(sys.argv) < 2 or not sys.argv[1] in (c[0] for c in commands): # or "-h" in sys.argv or "--help" in sys.argv:
        print_commands()
        sys.exit(1)

    cmd = next(c for c in commands if c[0] == sys.argv[1])
    cmd[2](sys.argv[2:])
    try:  # don't ask. but dont remove.
        sys.stdout.flush()
    except IOError:
        pass

if __name__ == "__main__":
    main()
