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


_driver_location = pjoin(os.path.abspath(os.path.curdir), "..",
                        "..", "control_scripts", "driver", "driver")


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


def _selftestrun(argdict):

    arglist = mm.argdict2list(argdict)

    outdir = pjoin("Validation", argdict['-ts'])

    mm.cleandir(outdir)

    # Run the command, put stdout + stderr into a file
    with open(pjoin(outdir, "stdout"), 'w') as hstdout:
        flag = subp.call(arglist + [ "-outdir", outdir],
                                stdout=hstdout,
                                stderr=subp.STDOUT)

    try:
        data = mm.parse_trace_file(pjoin(outdir, "trace"))
    except IOError:
        data = None

    return data, outdir, flag


def selftestrun(argdict):

    maxerrortol = 1e-3

    data, outdir, flag = _selftestrun(argdict)

    if flag != 0:
        fail_message(outdir)

    errors = data['error_norms']
    assert all(errors >= 0)
    maxerror = max(errors)

    if maxerror > maxerrortol:
        fail_message(outdir, maxerror)
        return False

    else:
        pass_message(outdir)
        return True


def main():

    argdicts = {
        "-binary" : [_driver_location],
        "-driver" : ["ll"],
        "-dt": [0.05],
        "-scale": [5],
        "-mesh": ["sq_square"],
        "-disable-ms" : [True],
        "-ts" : ["rk2", "rk4"],
        }

    passes = mm.parallel_parameter_sweep(selftestrun, argdicts,
                                         serial_mode=True)

    # argdicts = {
    #     "-binary" : [_driver_location],
    #     "-driver" : ["ll"],
    #     "-dt": [0.05],
    #     "-scale": [5],
    #     "-mesh": ["sq_square"],
    #     "-disable-ms" : [True],
    #     "-ts" : ["rk2", "rk4"],
    #     }

    # passes = mm.parallel_parameter_sweep(selftestrun, argdicts,
    #                                      serial_mode=True)

    # ??ds compare with implict!

    if all(passes):
        return 0
    else:
        return 1



if __name__ == "__main__":
    sys.exit(main())
