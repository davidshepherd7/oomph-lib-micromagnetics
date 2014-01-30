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

import oomphpy
import oomphpy.micromagnetics as mm


def pass_message(outdirname):
    print("PASSED", os.path.abspath(outdirname))


def fail_message(outdirname, maxerror=None, maxlteerror=None):

    absoutdirname = os.path.abspath(outdirname)

    if maxerror is not None:
        print("FAILED", absoutdirname, "max error of", maxerror)

    elif maxlteerror is not None:
        print("FAILED", absoutdirname, "max error in lte estimate",
              maxlteerror)

    else:
        print("FAILED", absoutdirname)
        with open(pjoin(outdirname, "stdout"), 'r') as hstdout:
            print(hstdout.read())


def constant_dt_test(exact, timestepper):
    maxerrortol = 0.1

    outdir = pjoin("Validation", exact + "_" + timestepper)

    mm.cleandir(outdir)

    arglist = ["ode",
               "-disable-explicit-solver-optimisations",
               "-outdir", outdir,
                "-dt", "0.05",
                 "-tmax", "4",
                  "-ts", timestepper,
                   "-exact", exact]

    # Run the command, put stdout + stderr into a file, put exact command
    # into another file
    flag = mm.run_driver(arglist, outdir)

    if flag != 0:
        fail_message(outdir)
        return False

    data = mm.parse_trace_file(pjoin(outdir, "trace"))
    maxerror = max(data['error_norms'])
    assert all(data['error_norms'] >= 0), exact +" "+ timestepper + str(data['error_norms'])


    if maxerror > maxerrortol:
        fail_message(outdir, maxerror)
        return False

    else:
        pass_message(outdir)
        return True


def adaptive_midpoint_test(exact, timestepper, predictor):

    tol = 1e-3

    outdir = pjoin("Validation", exact + "_adaptive" + timestepper
                   + "_" + predictor)

    mm.cleandir(outdir)

    arglist = ["ode",
               "-disable-explicit-solver-optimisations", # bugs? :(
               "-outdir", outdir,
               "-tol", str(tol),
               "-tmax", "4",
               "-ts", timestepper,
               "-mp-pred", predictor,
               "-exact", exact]
    flag = mm.run_driver(arglist, outdir)


    if flag != 0:
        fail_message(outdir)
        return False

    data = mm.parse_trace_file(pjoin(outdir, "trace"))
    maxerror = max(data['error_norms'])
    assert all(data['error_norms'] >= 0)

    # Bound the error with something proportional to the max value if it
    # gets large.
    maxvalue = max(it.chain(data['trace_values'], [1.0]))
    maxerrortol = 10 * tol * maxvalue


    # For second order polynomial we know lte is zero, check we got this
    # correct
    ltefail = False
    if exact == "poly2":
        maxlteerror = max(map(lambda a,b: a-b, data['LTE_norms'], it.repeat(0)))
        ltefail = (maxlteerror > 10e-8)

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

    passes = []

    exacts = ["sin", "cos", "poly3", "poly2"]
    timesteppers = ["rk2", "rk4", "midpoint-bdf", "bdf2"]

    for exact, timestepper in it.product(exacts, timesteppers):
        passes.append(constant_dt_test(exact, timestepper))


    mp_timesteppers = ["midpoint-bdf"]
    mp_preds = ["ebdf3", "rk2", "rk4"]

    for exact, timestepper, predictor in it.product(exacts, mp_timesteppers, mp_preds):
        passes.append(adaptive_midpoint_test(exact, timestepper, predictor))


    if not all(passes):
        return 1
    else:
        return 0


if __name__ == "__main__":
    sys.exit(main())
