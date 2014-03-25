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

# def ndt(dataset):
#     return len(dataset['dts'])


# def ndt_reldiff(ds):
#     if len(ds) != 2:
#         print("Not enough values")
#         return None

#     ndt1 = ndt(ds[0])
#     ndt2 = ndt(ds[1])

#     return abs(ndt1 - ndt2) / min(ndt1, ndt2)

def main():

    # What to run
    argdicts = {
        "-driver" : "ode",
        "-exact" : ["sin", "poly3", "poly2", "stiff_test"],
        "-ts" : ["bdf2", "midpoint-bdf", "tr"],
        "-tmax" : 5,
        "-tol" : 1e-5,
        "-disable-mm-opt" : True,
        }

    # ??ds not sure why we need "-disable-mm-opt",
    # but there's something stupid going on with mass matrix storage

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

    # Split into data for each time stepper (we assume that the order of
    # data is preserved here so that bdf_data[n] is the same exact solution
    # as imr_data[n] etc.).
    bdf_data = [d for d in datasets if d['-ts'] == 'bdf2']
    imr_data = [d for d in datasets if d['-ts'] == 'midpoint-bdf']
    tr_data = [d for d in datasets if d['-ts'] == 'tr']

    # Use bdf's nsteps as a maximum, tr and imr are more accurate than bdf2
    # so this should be true (unless bdf2's numerical damping has kicked in
    # and caused it to jump to a steady state too soon),
    max_steps = [1.2*len(d['times']) for d in bdf_data]

    # Check all the runs:
    imr_nsteps_ok = [tests.check_ndt_less_than(d, m)
                     for d, m in zip(imr_data, max_steps)]
    test_results.append(all(imr_nsteps_ok))

    tr_nsteps_ok = [tests.check_ndt_less_than(d, m)
                     for d, m in zip(tr_data, max_steps)]
    test_results.append(all(tr_nsteps_ok))


    if all(test_results):
        return 0
    else:
        return 1


if __name__ == "__main__":
    sys.exit(main())
