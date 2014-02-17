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

def ndt(dataset):
    return len(dataset['dts'])


def ndt_reldiff(ds):
    if len(ds) != 2:
        print("Not enough values")
        return None

    ndt1 = ndt(ds[0])
    ndt2 = ndt(ds[1])

    return abs(ndt1 - ndt2) / min(ndt1, ndt2)

def main():

    # What to run
    argdicts = {
        "-driver" : "ode",
        "-exact" : ["sin", "poly3", "poly2", "stiff_test"],
        "-ts" : ["bdf2", "midpoint-bdf"],
        "-tmax" : 5,
        "-tol" : 1e-5,
        # "-disable-explicit-solver-optimisations" : [True, False],
        }

    # Where it's going to end up
    base_outdir = os.path.abspath(pjoin(os.path.dirname(__file__), "Validation"))

    # Run
    err_codes, outdirs = mm.run_sweep(argdicts, base_outdir)

    # Get data
    datasets = list(filter(lambda d: d is not None, map(mm.parse_run, outdirs)))


    # Check things ran
    test_results = []
    test_results.append(all([e == 0 for e in err_codes]))

    # Check errors are small
    test_results.append(all([tests.check_error_norm(d, 0.05) for d in datasets]))

    # Check ndt are close for mp vs bdf2
    split_datasets = mm.split_up_stuff(datasets, ['-exact'])
    # map(lambda ds: assert(len(ds) == 2), split_datasets)
    ndt_rel_diffs = [ndt_reldiff(ds) for ds in split_datasets]
    test_results.append(all([ndt < 0.1 for ndt in ndt_rel_diffs if ndt is not None]))


    if all(test_results):
        return 0
    else:
        return 1


if __name__ == "__main__":
    sys.exit(main())
