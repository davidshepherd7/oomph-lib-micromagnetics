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
import itertools as it
import scipy as sp
import functools as ft
import glob

from os.path import join as pjoin
from functools import partial as par

# Make sure *this* versions oomphpy is in the path (before any other
# versions in other places)
sys.path.insert(1, pjoin(os.path.dirname(__file__), "../../etc"))
import oomphpy
import oomphpy.micromagnetics as mm

import fpdiff

# Old implementations of things:

def parallel_parameter_sweep(function, parameter_dictionary, serial_mode=False,
                             **kwargs):
    """Run function with all combinations of parameters in parallel using
    all available cores.

    parameter_lists should be a list of lists of parameters,
    """

    import multiprocessing

    # Generate a complete set of combinations of parameters
    parameter_sets = [dict(zip(parameter_dictionary, x))
                      for x in it.product(*parameter_dictionary.values())]

    # For debugging we often need to run in serial (to get useful stack
    # traces).
    if serial_mode:
        results = map(function, parameter_sets)

    else:
        # Run in all parameter sets in parallel
        results = multiprocessing.Pool(**kwargs).map(function, parameter_sets, 1)

    return results


def pass_message(outdirname):
    print("PASSED", os.path.abspath(outdirname))


def fail_message(outdirname, maxerror=None):

    absoutdirname = os.path.abspath(outdirname)

    if maxerror is not None:
        print("FAILED", absoutdirname, "max error of", maxerror)

    else:
        print("FAILED", absoutdirname)
        with open(pjoin(outdirname, "stdout"), 'r') as hstdout:
            print(hstdout.read())


def selftestrun(argdict):

    outdir = argdict.get("-outdir")
    if outdir is None:
        outdir = pjoin("Validation", argdict['-ts'] + str(argdict.get('-ms-method')))
        argdict['-outdir'] = outdir

    arglist, _, _ = mm.argdict2list(argdict)

    mm.cleandir(outdir)
    flag = mm.run_driver(arglist, outdir)

    try:
        data = mm.parse_trace_file(pjoin(outdir, "trace"))
    except IOError:
        data = None

    return data, outdir, flag


def test_error_norms(argdict):

    maxerrortol = 3e-3

    data, outdir, flag = selftestrun(argdict)

    if flag != 0:
        fail_message(outdir)
        return False

    errors = data['error_norms']
    assert all(errors >= 0)
    maxerror = max(errors)

    if maxerror > maxerrortol:
        fail_message(outdir, maxerror)
        return False

    else:
        pass_message(outdir)
        return True


def test_vs_implicit(argdict):

    implicit_dir = "validata"
    implicit_data = mm.parse_trace_file(pjoin(implicit_dir, "trace"))
    maxerrortol = 1e-3

    data, outdir, flag = selftestrun(argdict)

    if flag != 0:
        fail_message(outdir)
        return False

    mx_diffs = implicit_data['mean_mxs'] - data['mean_mxs']
    my_diffs = implicit_data['mean_mys'] - data['mean_mys']
    mz_diffs = implicit_data['mean_mzs'] - data['mean_mzs']

    maxerror = max(it.chain(mx_diffs, my_diffs, mz_diffs))

    if maxerror > maxerrortol:
        fail_message(outdir, maxerror)
        return False
    else:
        pass_message(outdir)
        return True


def main():

    # Parse arguments
    # ============================================================
    parser = argparse.ArgumentParser(description=main.__doc__,
    # Don't mess up my formating in the help message
    formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--generate-validata', action = "store_true")
    args = parser.parse_args()


    # Run tests
    # ============================================================

    # Test without magnetostatics by comparison with Mallinson solution
    noms_argdicts = {
        "-driver" : ["ll"],
        "-dt": [0.05],
        "-scale": [5],
        "-mesh": ["sq_square"],
        "-ref" : [2],
        "-ms-method" : ["disabled"],
        "-ts" : ["rk2", "rk4"],
        "-tmax" : [3],
        }
    noms_passes = parallel_parameter_sweep(test_error_norms, noms_argdicts,
                                           serial_mode=True)

    # If requested then run fully implicit version to generate validata for
    # magnetostatics test.
    if args.generate_validata:
        implicit_argdict = {
            "-driver" : "ll",
            "-dt": 0.01,
            "-scale": 5,
            "-ref" : [2],
            "-mesh": "sq_square",
            "-ms-method" : ["decoupled"],
            "-ts" : "midpoint-bdf",
            "-fd-jac" : 1,
            "-outdir" : "validata",
            }
        selftestrun(implicit_argdict)


    # Test with magnetostatics by comparison with (fully) implicit
    # timesteppers. Implicit needs to work the same as decoupled so that
    # adaptive midpoint works.
    ms_argdicts = {
        "-driver" : ["ll"],
        "-dt": [0.01],
        "-scale": [5],
        "-ref" : [2],
        "-mesh": ["sq_square"],
        "-ms-method" : ["decoupled", "implicit"],
        "-ts" : ["rk2"], # Don't run rk4 because 4th order makes it
                         # hard to compare with implicit methods (no
                         # A-stable 4th order methods).
        }

    # Now run explicit ones and compare
    ms_passes = parallel_parameter_sweep(test_vs_implicit, ms_argdicts,
                                         serial_mode=True)

    if all(noms_passes) and all(ms_passes):
        return 0
    else:
        return 1



if __name__ == "__main__":
    sys.exit(main())
