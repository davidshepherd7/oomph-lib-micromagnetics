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


def try_run(args, dt, dir_naming_args, root_outdir,
             maxallowederror = 0.01, maxallowedangle = sp.pi/4):

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

    success = (err_code == 0 and maxerror < maxallowederror
                and maxangle < maxallowedangle)

    if success:
        mm.okprint("Succedded in", data['-outdir'])
    else:
        mm.badprint(pjoin(data['-outdir'], "stdout:1:1:"), "FAILED")

    return success, data




def locate_stable_point(args, refine, dt_guess, **kwargs):

    # Catch failures immediately
    if dt_guess is None:
        return None

    # Otherwise start from the initial guess
    else:
        dt = dt_guess


    success = False

    args.update({'-ref' : refine})

    # If it gets below this dt then give up
    while not success and (dt >= 1e-6):

        # Try a run
        success, data = try_run(args, dt, **kwargs)

        # On fail half the step and try again
        if not success:
            dt = dt/2


    if not success:
        # if we get here then it totally failed
        return None
    else:
        return dt


def locate_stable_points(args, refines_list, dt_guess, dir_naming_args,
                         root_outdir, nbisection=0):

    # Also name by refinement and dt
    results_naming_args = copy.copy(dir_naming_args)
    dir_naming_args = dir_naming_args + ['-ref', '-dt']

    results_naming_values = [str(args[k]) for k in sorted(results_naming_args)]
    results_name = pjoin(root_outdir, 'stable_dts_'+'_'.join(results_naming_values))

    print("writing results to", results_name)

    # # Run for easiest case
    # dt = locate_stable_point(args, refine=refines_list[0],
    #                          dt_guess=dt_guess, dir_naming_args=dir_naming_args,
    #                          root_outdir=root_outdir)

    # # Write to file
    # with open(results_name, 'a') as result_file:
    #     result_file.write(str(refines_list[0]) + " " + str(dt) + "\n")

    dt_ranges = []
    dt = dt_guess

    # And the rest
    for ref in refines_list:

        # find stable point for this refine, use previous dt as initial guess
        dt = locate_stable_point(args,
                                 refine=ref,
                                 dt_guess=dt,
                                 dir_naming_args=dir_naming_args,
                                 root_outdir=root_outdir)

        # Refine for more accurate dt value using binary search
        if dt < dt_guess:
            dt_range = refine_stable_point(args, dt, 2*dt, nbisection,
                                        dir_naming_args=dir_naming_args,
                                        root_outdir=root_outdir)
        else:
            dt_range = (dt_guess, dt_guess)

        dt_ranges.append(dt_range)

        # Output
        print("Final range", dt_range)
        with open(results_name, 'a') as result_file:
            result_file.write(str(ref)
                              + " " + str(dt_range[0])
                              + " " + str(dt_range[1])
                              + "\n")

        # Next refinement should use the maximum possibly stable value as
        # its initial dt guess
        _, dt = dt_range




    return dt_ranges


def refine_stable_point(args, dt_min, dt_max, nrefines, **kwargs):

    # Recursion base case
    if nrefines == 0:
        return dt_min, dt_max

    print("searching in", dt_min, dt_max, "for", nrefines, "times")

    # See if the middle dt is stable
    dt_mid = (dt_min + dt_max)/2
    success, data = try_run(args, dt_mid, **kwargs)

    # Recurse with appropriate dt min/max and one less nrefine
    if success:
        return refine_stable_point(args, dt_mid, dt_max, nrefines-1, **kwargs)
    else:
        return refine_stable_point(args, dt_min, dt_mid, nrefines-1, **kwargs)


def main():

    parser = argparse.ArgumentParser(description=main.__doc__,

    # Don't mess up my formating in the help message
    formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--clean', action='store_true', help='delete old data')
    parser.add_argument('--serial', action='store_true',
                         help="Don't do parallel sweeps")
    parser.add_argument('--outdir', help='delete old data')
    parser.add_argument('--quick', action='store_true',
                        help='quick version of parameters for testing')
    args = parser.parse_args()

    if args.outdir is None:
        args.outdir = os.path.abspath('../experiments/intermag_rk4')

    argsdict = {
                '-driver' : ['ll', 'llg'],
                '-tmax' : 4,
                '-hlib-bem' : 0,
                '-renormalise' : 0,
                '-mesh' : ['ut_sphere'], # 'sq_square'],
                '-ms-method' : ['decoupled'],
                '-solver' : 'som-gmres',
                '-prec' : 'som-main-ilu-1',
                '-ts' : ['rk4'],
                '-scale' : 2,
                '-fd-jac' : True,
                '-damping' : [0.1],
                '-check-angles' : 1,
                }
    refines = [1, 2, 3, 4]

    # Only do a couple of things for a test run
    if args.quick:
        argsdict['-ms-method'] = 'disabled'
        argsdict['-damping'] = 1.0
        argsdict['-ts'] = 'euler'
        refines = [1, 2, 3]


    # The exact function to run
    f = par(locate_stable_points,
            refines_list = refines,
            dt_guess=0.1,
            dir_naming_args=mm.argdict_varying_args(argsdict),
            root_outdir=args.outdir,
            nbisection=2)

    # Create list of all possible combinations of args
    argsets = mm.product_of_argdict(argsdict)

    # Filter out impossible combinations, bit hacky...
    def is_bad_combination(data):
        if data['-ts'] == 'rk2' or data['-ts'] == 'euler'  or data['-ts'] == 'rk4':
            if data['-driver'] == 'llg' or data['-ms-method'] == 'implicit':
                return True
        elif data['-ts'] == 'midpoint-bdf' or data['-ts'] == 'bdf2':
            if data['-driver'] == 'll':
                return True

        return False

    argsets = [a for a in argsets if not is_bad_combination(a)]


    dts = list(mm.parallel_map(f, argsets, serial_mode=args.serial))

    print(dts)

    return 0

# If this script is run from a shell then run main() and return the result.
if __name__ == "__main__":
    sys.exit(main())
