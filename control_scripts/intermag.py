#!/usr/bin/env python3

# Future proofing
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

# Imports from main libraries
import subprocess as subp
from multiprocessing import Pool
import multiprocessing
import sys
import argparse
import os
import os.path
import hashlib
import shutil
import re

import itertools as it
import functools as ft
import scipy as sp
import matplotlib.pyplot as plt


# Imports for specific functions
from functools import partial as par
from os.path import join as pjoin
from glob import glob


import parse


library_folder = os.path.abspath("../")
driver_folder = pjoin(library_folder, 'control_scripts', 'driver')
root_outdir = '../experiments/intermag'

def locate_stable_point(refine, dt_guess, other_args):

    maxallowederror = 0.1
    dt = dt_guess
    success = False

    while (dt >= 1e-6):

        print("\nTrying with dt =", dt, "ref =", refine)

        outdir = pjoin(root_outdir, str(refine) + '-' + str(dt))
        try:
            os.makedirs(outdir)
        except OSError:
            pass

        # Try to run it
        with open(pjoin(outdir, "stdout"), 'w') as stdout_file:

            status = subp.call([pjoin(driver_folder, 'driver'),
                               'll',
                                '-disable-ms',
                                 '-tmax', '2',
                                 '-error-norm-limit', str(maxallowederror),
                                 '-outdir', outdir,
                                 '-ref', str(refine),
                                 '-dt', str(dt)]
                                 + other_args,
                                 stdout=stdout_file,
                                 stderr=subp.STDOUT)


        errs = parse.parse_trace_file(pjoin(outdir, 'trace'))['error_norms']

        assert errs[0] >= 0, errs

        maxerror = max(errs)

        # If it didn't crash then check output
        if status == 0:

            if maxerror < maxallowederror:
                 print("Succedded with max error =", maxerror, "dt =", dt)
                 return dt
            else:
                print("Failed due to max error =", maxerror, ">", maxallowederror)
                dt = dt/2

        else:
            print("Failed due to crash with max error =", maxerror)
            dt = dt/2


def main():

    parser = argparse.ArgumentParser(description=main.__doc__,

    # Don't mess up my formating in the help message
    formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--clean', action='store_true', help='delete old data')
    parser.add_argument('--outdir', '-o', help='Output dir')
    parser.add_argument('--mesh', '-m', help='Choose mesh type')
    parser.add_argument('--implicit', '-i', action='store_true',
                         help='Use implicit time integrator')

    args = parser.parse_args()

    if args.outdir is not None:
        root_outdir = args.outdir

    if args.clean:
        shutil.rmtree(root_outdir)

    if args.mesh == "sphere":
        mesh_args = ['-mesh', 'ut_sphere', '-scale', '2.5']
    elif args.mesh == "square":
        mesh_args = ['-mesh', 'sq_square', '-scale', '5']

    else:
        sys.stderr.write("Unrecognised mesh" + str(args.mesh))
        exit(2)


    if args.implicit:
        ts_args = ['-ts', 'midpoint-bdf', '-fd-jac']
    else:
        ts_args = ['-ts', 'rk2']



    # Make sure micromag library is up to date
    library_folder = os.path.abspath("../")
    print("Building and installing libraries from", library_folder)
    subp.check_call(['make', '--silent', '--keep-going',
                     'LIBTOOLFLAGS=--silent'], cwd=library_folder)
    subp.check_call(['make', 'install', '--silent', '--keep-going',
                     'LIBTOOLFLAGS=--silent'], cwd=library_folder)

    # Make sure the driver binaries are up to date
    subp.check_call(['make', '--silent', '--keep-going',
                    'LIBTOOLFLAGS=--silent'], cwd=driver_folder)

    # Run
    refs = [1, 2, 3, 4, 5]
    dts = []
    dt_found = 1
    for ref in refs:
        dt_found = locate_stable_point(ref, dt_found, mesh_args + ts_args)
        dts.append(dt_found)

    # Write to file
    with open(pjoin(root_outdir, 'results'), 'w') as result_file:
        for r, dt in zip(refs, dts):
            result_file.write(str(r) + " " + str(dt))

    plt.plot(refs, dts)
    plt.xlabel('refinement')
    plt.ylabel('max stable dt found')
    plt.show()

    return 0

# If this script is run from a shell then run main() and return the result.
if __name__ == "__main__":
    sys.exit(main())
