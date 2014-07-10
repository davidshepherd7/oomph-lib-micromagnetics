#!/usr/bin/env python3

# Python 2/3 compatability
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import sys
import argparse
import os
import os.path

from os.path import join as pjoin

# Make sure *this* versions oomphpy is in the path (before any other
# versions in other places)
sys.path.insert(1, pjoin(os.path.dirname(__file__), "../../etc"))
import oomphpy
import oomphpy.micromagnetics as mm
import oomphpy.tests as tests


def main():

    # Look for parallel in args
    parser = argparse.ArgumentParser()
    parser.add_argument('--parallel', action = "store_true")
    args = parser.parse_args()

    # What to run
    argdicts = {
        "-driver" : 'llgode',
         '-exact' : "ll",
         '-ts' : ["rk2", "midpoint-bdf"],
         '-dt' : 0.01,
         '-damping' : 0.5,
         '-h-app' : 'minus_z',
        }

    # Where it's going to end up
    base_outdir = os.path.abspath(pjoin(os.path.dirname(__file__), "Validation"))

    # Run
    err_codes, outdirs = mm.run_sweep(argdicts, base_outdir,
                                      parallel_sweep=args.parallel)

    # Get data
    datasets = list(map(mm.parse_run, outdirs))

    # Check all errors are small
    ok = all([tests.check_error_norm(d, 1e-4) for d in datasets])

    ran = all([e == 0 for e in err_codes])

    if ran and ok:
        return 0
    else:
        return 1


if __name__ == "__main__":
    sys.exit(main())
