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

# Make sure *this* versions oomphpy is in the path (before any other
# versions in other places)
sys.path.insert(1, pjoin(os.path.dirname(__file__), "../etc"))
import oomphpy
import oomphpy.micromagnetics as mm



def locate_stable_point(refine, dt_guess, other_args, root_outdir):

    maxallowederror = 0.1
    dt = dt_guess
    success = False

    while (dt >= 1e-6):

        print("\nTrying with dt =", dt, "ref =", refine)

        # Main list of args
        args = {'-driver' : 'll',
                '-tmax' : 2,
                '-error-norm-limit' : maxallowederror,
                '-ref' : refine,
                '-dt' : dt,
                '-hlib-bem' : 0,
                }

        # Merge in args given
        args.update(other_args)

        # Args to use when generating a dir name
        dir_naming_args = ['-ref', '-dt']

        # Try to run it
        err_code, outdir = mm._run(args, root_outdir, dir_naming_args)

        # If it didn't crash then check output
        if err_code == 0:

            # Parse output files
            data = mm.parse_run(outdir)
            errs = data['error_norms']
            assert errs[0] >= 0

            maxerror = max(errs)

            if maxerror < maxallowederror:
                 print("Succedded with max error =", maxerror, "dt =", dt)
                 return dt
            else:
                print("Failed due to max error =", maxerror, ">", maxallowederror)
                dt = dt/2

        else:
            print("Failed due to crash")
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
    parser.add_argument('--use-hms', action='store_true', help='Include magnetostatics')


    args = parser.parse_args()

    if args.outdir is not None:
        root_root_outdir = args.outdir
    else:
        root_root_outdir = os.path.abspath('../experiments/intermag')

    if args.clean:
        shutil.rmtree(root_outdir)


    additional_args = {}

    if args.mesh == "sphere":
        additional_args.update({'-mesh' : 'ut_sphere', '-scale' : '2.5'})
    elif args.mesh == "square":
        additional_args.update({'-mesh' : 'sq_square', '-scale' : '5'})
    else:
        sys.stderr.write("Unrecognised mesh " + str(args.mesh))
        exit(2)


    if args.implicit:
        additional_args.update({'-ts' : 'midpoint-bdf', '-fd-jac' : True})
    else:
        additional_args.update({'-ts': 'rk2'})


    if args.use_hms:
        additional_args.update({'-solver' : 'som-gmres', '-prec' : 'som-main-exact'})
    else:
        additional_args.update({'-disable-ms' : True})


    # Construct root dir
    args_dirname = '{}_impl{}_hms{}'.format(args.mesh, args.implicit, args.use_hms)
    root_outdir = pjoin(root_root_outdir, args_dirname)

    # Run
    refs = [1, 2, 3, 4, 5]
    dts = []
    dt_found = 1
    for ref in refs:
        dt_found = locate_stable_point(ref, dt_found, additional_args,
                                       root_outdir)
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
