#!/usr/bin/env python3

# Python 2/3 compatability
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import sys
import argparse

import os
import shutil
import os.path
import scipy as sp
import itertools as it

from os.path import join as pjoin

# Make sure *this* versions oomphpy is in the path (before any other
# versions in other places)
sys.path.insert(1, pjoin(os.path.dirname(__file__), "../../etc"))
import oomphpy
import oomphpy.micromagnetics as mm
import oomphpy.tests as tests


def adaptive_midpoint_test(exact, timestepper, predictor):


    # ??ds still not testing ltes for other odes... probably should


    if maxerror > maxerrortol:
        fail_message(outdir, maxerror)
        return False

    elif ltefail:
        fail_message(outdir, maxlteerror=maxlteerror)
        return False

    else:
        pass_message(outdir)
        return True


def main():

    # Look for parallel in args
    parser = argparse.ArgumentParser()
    parser.add_argument('--parallel', action = "store_true")
    args = parser.parse_args()

    # Where it's going to end up
    base_outdir = os.path.abspath(pjoin(os.path.dirname(__file__), "Validation"))

    # constant dt
    # ============================================================

    argdicts_const_dt = {
        '-driver': "ode",
        "-disable-mm-opt" : True,
        "-dt": "0.05",
        "-tmax": "4",
        "-ts": ["rk2", "rk4", "midpoint-bdf", "bdf2", "tr"],
        "-exact": ["sin", "cos", "poly3", "poly2"],
        }

    # Run const
    err_codes_const_dt, outdirs_const_dt = \
      mm.run_sweep(argdicts_const_dt, pjoin(base_outdir, "const_dt"),
                   parallel_sweep=args.parallel)

    # Check they all ran without crashing
    const_ran = all(e == 0 for e in err_codes_const_dt)

    # Parse output and check error norms
    const_tests = [tests.check_error_norm(data, tol=0.1)
                   for data in map(mm.parse_run, outdirs_const_dt)]


    # varying dt
    # ============================================================

    argdicts_var_dt = {
        '-driver': "ode",
        "-disable-mm-opt" : True,
        "-dt-initial": 1e-6,
        "-tol": 1e-3,
        "-tmax": "4",
        "-ts": ["midpoint-bdf"],
        "-exact": ["sin", "cos", "poly3", "poly2"],
        "-mp-pred" : ["ebdf3", "rk2", "rk4"],
        }

    # Run var
    err_codes_var_dt, outdirs_var_dt = \
      mm.run_sweep(argdicts_var_dt, pjoin(base_outdir, "var_dt"),
                   parallel_sweep=args.parallel)

    # Check they all ran without crashing
    var_ran = all(e == 0 for e in err_codes_var_dt)

    # Parse output
    datasets = list(map(mm.parse_run, outdirs_var_dt))

    # check error norms
    var_tests = [tests.check_error_norm_relative(data, tol=1e-2)
                 for data in datasets]

    # For second order polynomial we know lte is zero, check we got these cases
    # correct.
    poly2_data = [d for d in datasets if d['exact'] == 'poly2']
    lte_test = [tests.check_error_norm(d, tol=1e-7) for d in poly2_data]


    # return
    # ============================================================

    if const_ran and all(const_tests) \
      and var_ran and all(var_tests) and all(lte_test):
        return 0
    else:
        return 1


if __name__ == "__main__":
    sys.exit(main())
