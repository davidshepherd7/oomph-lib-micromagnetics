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
        "-ts" : ["bdf2", "midpoint-bdf", "tr", "bdf1"],
        "-tmax" : 10,
        "-tol" : 1e-4,
        "-disable-mm-opt" : True,
        "-always-write-trace" : 1, # Otherwise we get wrong ndts by counting len(dts)
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


    # Use bdf2's nsteps as a maximum, tr and imr are more accurate than
    # bdf2 so this should be true (unless bdf2's numerical damping has
    # kicked in and caused it to jump to a steady state too soon), (we
    # assume that the order of data is preserved here so that bdf_data[n]
    # is the same exact solution as imr_data[n] etc.).
    bdf2_data = [d for d in datasets if d['-ts'] == 'bdf2']



    for ts in argdicts['-ts']:

        # bdf1 sucks (first order) so it needs far more steps, do it
        # manually.
        if ts == "bdf1":
            max_err = 0.4
            max_steps = [550, 3800, 1050, 70]
        else:
            max_err = 0.07
            max_steps = [1.3*len(d['times']) for d in bdf2_data]

        ts_data = [d for d in datasets if d['-ts'] == ts]

        # Check errors are small
        test_results.append(all([tests.check_error_norm(d, max_err) for d in ts_data]))

        # Check all the runs:
        nsteps_ok = [tests.check_ndt_less_than(d, m)
                     for d, m in zip(ts_data, max_steps)]
        test_results.append(all(nsteps_ok))


    if all(test_results):
        return 0
    else:
        return 1


if __name__ == "__main__":
    sys.exit(main())
