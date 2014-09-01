#!/usr/bin/env python3

# Python 2/3 compatability
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import sys
import argparse
import os
import os.path

import itertools as it

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
        "-driver" : ["llg"],
        "-ts" : ["midpoint-bdf", "bdf2"],
        "-ms-method" : ["implicit", "decoupled", "disabled"],
        "-solver" : "gmres",
        "-matrix-type" : "som",
        "-prec" : "som-main-exact",
        "-tmax" : 100,
        "-tol" : 0.01,
        "-mp-pred" : "ebdf3",
        # "-initial-m" : "exactly_z",
        }

    # Where it's going to end up
    base_outdir = os.path.abspath(pjoin(os.path.dirname(__file__), "Validation"))

    # Run
    err_codes, outdirs = mm.run_sweep(argdicts, base_outdir,
                                      serial_mode=not args.parallel)

    # Get data
    datasets = list(map(mm.parse_run, outdirs))

    # Check things
    ran = all((e == 0 for e in err_codes))
    t1 = all([tests.check_ndt_less_than(d, 80) for d in datasets
              if d['-ms-method'] != "decoupled"])
    t2 = all([tests.check_ndt_less_than(d, 115) for d in datasets
              if d['-ms-method'] == "decoupled"])

    # Decoupled ms introduces wiggles which put an effective cap on the step size.

    if all([ran, t1, t2]):
        return 0
    else:
        return 1


if __name__ == "__main__":
    sys.exit(main())
