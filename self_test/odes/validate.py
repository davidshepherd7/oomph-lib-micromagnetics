#!/usr/bin/env python

# Python 2/3 compatability
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import sys
import argparse

import os
import shutil
import os.path
import subprocess as subp
import itertools as it
import scipy as sp

from os.path import join as pjoin

import oomphpy
import oomphpy.micromagnetics as mm


driver = pjoin(os.path.abspath(os.path.curdir), "..",
                "..", "control_scripts", "driver", "driver")


def cleandir(dirname):
    """(Re)make a directory called dirname.
    """

    # If it exists then delete all files (won't touch subdirs though) and
    # the folder itself. This will fail if we have any subdirs (for safety).
    if os.path.isdir(dirname):
        for f in os.listdir(dirname):
            os.unlink(pjoin(dirname, f))
        os.rmdir(dirname)

    # Make the directory and any parents needed
    os.makedirs(dirname)


def constant_dt_test(exact, timestepper):
    maxerrortol = 0.1

    outdir = pjoin("Validation", exact + "_" + timestepper)

    cleandir(outdir)

    # Run the command, put stdout + stderr into a file
    with open(pjoin(outdir, "stdout"), 'w') as hstdout:
        flag = subp.call([driver, "ode",
                         "-disable-explicit-solver-optimisations", # bugs? :(
                         "-outdir", outdir,
                         "-dt", "0.05",
                         "-tmax", "4",
                         "-ts", timestepper,
                         "-exact", exact],
                        stdout=hstdout,
                        stderr=subp.STDOUT)


    data = mm.parse_trace_file(pjoin(outdir, "trace"))
    maxerror = max(data['error_norms'])

    # If it crashed then cat the stdout file
    if flag != 0:
        print("FAILED", exact, timestepper, "run did not finish")
        with open(pjoin(outdir, "stdout"), 'r') as hstdout:
            print(hstdout.read())

    elif maxerror > maxerrortol:
        print("FAILED", exact, timestepper, "max error of", maxerror, ">",
              maxerrortol)

    else:
        print("PASSED", exact, timestepper)


def adaptive_midpoint_test(exact, timestepper, predictor):

    tol = 1e-3
    maxerrortol = 10* tol

    outdir = pjoin("Validation", exact + "_adaptive" + timestepper
                   + "_" + predictor)

    cleandir(outdir)

    # Run the command, put stdout + stderr into a file
    with open(pjoin(outdir, "stdout"), 'w') as hstdout:
        flag = subp.call([driver, "ode",
                         "-disable-explicit-solver-optimisations", # bugs? :(
                         "-outdir", outdir,
                         "-tol", str(tol),
                         "-tmax", "4",
                         "-ts", timestepper,
                         "-mp-pred", predictor,
                         "-exact", exact],
                        stdout=hstdout,
                        stderr=subp.STDOUT)

    data = mm.parse_trace_file(pjoin(outdir, "trace"))
    maxerror = max(data['error_norms'])


    # If it crashed then cat the stdout file
    if flag != 0:
        print("FAILED", exact, "adaptive-" + timestepper, "pred-" + predictor)
        with open(pjoin(outdir, "stdout"), 'r') as hstdout:
            print(hstdout.read())

    elif maxerror > maxerrortol:
        print("FAILED", exact, timestepper, "max error of", maxerror, ">",
              maxerrortol)

    else:
        print("PASSED", exact, "adaptive-" + timestepper, "pred-" + predictor)


def main():

    exacts = ["sin", "cos", "poly3", "poly2"]
    timesteppers = ["rk2", "rk4", "midpoint", "midpoint-bdf", "bdf2"]


    for exact, timestepper in it.product(exacts, timesteppers):
        constant_dt_test(exact, timestepper)


    mp_timesteppers = ["midpoint", "midpoint-bdf"]
    mp_preds = ["ebdf3", "rk2", "rk4"]
    for exact, timestepper, predictor in it.product(exacts, mp_timesteppers, mp_preds):
        adaptive_midpoint_test(exact, timestepper, predictor)


    return 0


if __name__ == "__main__":
    sys.exit(main())
