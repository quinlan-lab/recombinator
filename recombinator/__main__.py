from __future__ import print_function, absolute_import
import warnings

warnings.simplefilter('ignore')
from . import recombinator
from . import cohort_plots
from . import __version__

import sys

commands = [
        ('recombinator', 'find informative sites and crossovers', recombinator.main),
        ('filter', 'call crossovers from informatives sites', recombinator.filter_main),
        ('cohort-plots', 'make plots summarizing crossovers across a cohort', cohort_plots.main),
        ]

def print_commands():
    print("recombinator version: %s" % __version__)
    print("""
commands:
%s""" % "\n".join(("%-14s: %s" % (c[0], c[1])) for c in commands))


def main():
    if len(sys.argv) < 2 or not sys.argv[1] in (c[0] for c in commands) or "-h" in sys.argv or "--help" in sys.argv:
        print_commands()
        sys.exit(1)

    cmd = next(c for c in commands if c[0] == sys.argv[1])
    cmd[2](sys.argv[2:])

if __name__ == "__main__":
    main()
