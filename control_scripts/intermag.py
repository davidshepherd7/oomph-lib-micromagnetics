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
import copy

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



def locate_stable_point(args, refine, dt_guess, dir_naming_args,
                        root_outdir):

    maxallowederror = 0.01
    maxallowedangle = sp.pi/4 # radians

    dt = dt_guess
    success = False

    args.update({'-ref' : refine})

    while (dt >= 1e-6):

        args.update({'-dt' : dt})

        # Try to run it
        err_code, outdir = mm._run(args, root_outdir, dir_naming_args,
                                    quiet=True)


        # Parse output files
        data = mm.parse_run(outdir)
        if data is not None:
            errs = data['m_length_error_means']
            assert errs[0] >= 0
            maxerror = max(errs)

            angles = data['max_angle_errors']
            maxangle = max(angles)

            maxtime = data['times'][-1]
            nsteps = len(data['times'])

        else:
            maxerror = sp.inf
            maxangle = sp.inf
            maxtime = 0
            nsteps = 0


        # If it didn't crash then check output
        if (err_code != 0
            or maxerror > maxallowederror
            or maxangle > maxallowedangle):
            mm.badprint(pjoin(outdir, "stdout:1:1:"), "FAILED",
                        maxerror, maxangle)

            dt = dt/2

        else:
            mm.okprint("Succedded in", os.path.relpath(outdir, root_outdir))
            return dt


def locate_stable_points(args, refines_list, dt_guess, dir_naming_args,
                         root_outdir):

    # Also name by refinement and dt
    results_naming_args = copy.copy(dir_naming_args)
    dir_naming_args = dir_naming_args + ['-ref', '-dt']

    results_naming_values = sorted([str(args[k]) for k in results_naming_args])
    results_name = pjoin(root_outdir, 'stable_dts_'+'_'.join(results_naming_values))

    print("writing results to", results_name)

    # Run for easiest case
    dt = locate_stable_point(args, refine=refines_list[0],
                             dt_guess=dt_guess, dir_naming_args=dir_naming_args,
                             root_outdir=root_outdir)

    # Write to file
    with open(results_name, 'a') as result_file:
        result_file.write(str(refines_list[0]) + " " + str(dt) + "\n")

    dts = [dt]

    # And the rest
    for ref in refines_list[1:]:

        # find stable point for this refine, use previous dt as initial guess
        dt = locate_stable_point(args,
                                 refine=ref,
                                 dt_guess=dts[-1],
                                 dir_naming_args=dir_naming_args,
                                 root_outdir=root_outdir)
        dts.append(dt)

        # Write to file
        with open(results_name, 'a') as result_file:
            result_file.write(str(ref) + " " + str(dt) + "\n")


    return dts


def main():

    parser = argparse.ArgumentParser(description=main.__doc__,

    # Don't mess up my formating in the help message
    formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--clean', action='store_true', help='delete old data')
    parser.add_argument('--serial', action='store_true',
                         help="Don't do parallel sweeps")
    args = parser.parse_args()

    root_root_outdir = os.path.abspath('../experiments/intermag')

    argsdict = {
                '-driver' : 'll',
                '-tmax' : 2,
                '-hlib-bem' : 0,
                '-renormalise' : 0,
                '-mesh' : ['ut_sphere'], # 'sq_square'],
                '-ms-method' : ['implicit', 'disabled', 'decoupled'],
                '-solver' : 'som-gmres',
                '-prec' : 'som-main-exact',
                '-ts' : ['midpoint-bdf', 'rk2'],
                '-scale' : 2,
                '-fd-jac' : True,
                }


    # The exact function to run
    f = par(locate_stable_points, refines_list = [1, 2, 3, 4, 5],
            dt_guess=0.1, dir_naming_args=mm.argdict_varying_args(argsdict),
            root_outdir=root_root_outdir)

    # Run for all combinations of args
    argsets = mm.product_of_argdict(argsdict)
    dts = list(mm.parallel_map(f, argsets, serial_mode=args.serial))

    print(dts)

    # # Write to file
    # with open(pjoin(root_outdir, 'results'), 'w') as result_file:
    #     for r, dt in zip(refs, dts):
    #         result_file.write(str(r) + " " + str(dt))

    # plt.plot(refs, dts)
    # plt.xlabel('refinement')
    # plt.ylabel('max stable dt found')
    # plt.show()

    return 0

# If this script is run from a shell then run main() and return the result.
if __name__ == "__main__":
    sys.exit(main())
