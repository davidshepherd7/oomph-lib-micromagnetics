#!/usr/bin/env python

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

# My code
import oomphpy
import oomphpy.micromagnetics as mm

# Imports for specific functions
from functools import partial as par
from os.path import join as pjoin
from glob import glob


# possible improvements:

# On error put stdout file somewhere useful?


greenColour = '\033[01;32m'
redColour = '\033[01;31m'
endColour = '\033[0m'


def execute_oomph_driver(args_dict, output_root):

    # Convert the dict to a list
    arglist, binary, mpi = mm.argdict2list(args_dict)

    # Construct an output directory name based on inputs if one has not
    # been specified.
    outdir = args_dict.get('-outdir')
    if outdir is not None:
        error("Don't specify an outdir, it will be automatically generated")

    # Create a hash of the inputs and use it to label the folder
    h = hashlib.sha224( ''.join([str(v) for _, v in args_dict.items()]).encode())
    outdir = pjoin(output_root,"results_" + h.hexdigest())
    arglist = arglist + ['-outdir', outdir]

    # Make sure the directory is empty and exists
    mm.cleandir(outdir)

    # Run and write stdout, stderr, command details
    err_code  = mm.run_driver(arglist, outdir, binary=binary, mpi_command=mpi)

    # Compress the output data files (but not info or trace).
    for datfilename in glob(pjoin(outdir, "*.dat")):
        subp.check_call(['gzip', datfilename])

    # Report result
    command = ' '.join(arglist)
    if err_code == 0:
        print(greenColour, command, endColour)
    else:
        # Print failure message to screen
        print('\n', redColour, command , endColour)
        print(redColour, "FAILED with exit code", err_code, "see",
              pjoin(outdir, "stdout"), endColour)

        # Print failure message into file in ouput directory
        print('This run failed!', file=open(pjoin(outdir, "FAILED"), 'w'))

    return 0


def standard_sweep(parameter_set, cleanup, serial_mode=False, no_build=False):

    if parameter_set == 'script_test':
        args_dict = {
            '-driver' : ["ode"],
            '-dt' : [1e-4],
            '-tmax' : [1.0],
            '-tol' : [1e-3],
            '-ref' : [2],
            }
    elif parameter_set == "compare-implicitness-implicit":
         args_dict = {
            '-binary' : ["./llg_driver/llg_driver"],
            '-tmax' : [20],
            '-ts' : ["bdf2"],
            '-mesh' : ['many_ut_square'],
            '-ref' : [3, 4],
            '-tol' : [1e-1, 1e-3, 1e-5],
            '-dt' : [1e-6],
            '-implicit-ms' : [True],
            '-solver' : ['som-gmres'],
            '-prec' : ['som-main-exact']
            }

    elif parameter_set == "check-single-ele-mesh":
         args_dict = {
            '-binary' : ["./llg_driver/llg_driver"],
            '-tmax' : [10],
            '-ts' : ["bdf2", "rk4"],
            '-mesh' : ['single-element'],
            '-dt' : [1e-2, 1e-1, 1e-3],
            '-resi' : ['ll'],
            '-fd-jac' : [True],
            }


    elif parameter_set == "coarse-blocked-ut-preconditioner":
         args_dict = {
            '-binary' : ["./llg_driver/llg_driver"],
            '-tmax' : [10.0],
            '-ts' : ["bdf2"],
            '-mesh' : ['sq_square', 'ut_sphere'],
            '-dt' : [0.5, 0.1],
            '-ref' : [2, 4, 5],
            '-implicit-ms' : [True],
            '-solver' : ['som-gmres'],
            '-prec' : ['som-main-blockut'],
            '-blocking' : ['group-m-phi-phi-boundary'],
            '-scale' : [1, 1000],
            }

    elif parameter_set == "blah":
         args_dict = {
            '-binary' : ["./semi_implicit_mm_driver/semi_implicit_mm_driver"],
            '-tmax' : [20],
            '-ts' : ["rk2"],
            '-dt' : [0.1, 1e-2, 1e-3, 1e-4],
            '-scale' : [10],
            '-resi' : ['ll'],
            }


    elif parameter_set == "blah2":
         args_dict = {
            '-binary' : ["./llg_driver/llg_driver"],
            '-tmax' : [20],
            '-ts' : ["bdf2"],
            '-tol' : [1e-3, 1e-4],
            '-scale' : [10],
            '-implicit-ms' : [True],
            '-solver' : ['som-gmres'],
            '-prec' : ['som-main-blockut'],
            }

    elif parameter_set == "implicit-vs-explicit-ms":
         args_dict = {
            '-binary' : ["./driver/driver"],
            '-driver' : ['llg'],
            '-tmax' : [10.0],
            '-ts' : ["bdf2"],
            '-mesh' : ['sq_cubeoid'],
            '-tol' : [1e-3],
            '-ref' : [0, 1],
            '-decoupled-ms' : [True, False],
            '-solver' : ['som-gmres'],
            '-prec' : ['som-main-blockut'],
            '-blocking' : ['group-m-phi-phi-boundary'],
            }

    elif parameter_set == "decoupled-ms-debug":
         args_dict = {
            '-binary' : ["./driver/driver"],
            '-driver' : ['llg'],
            '-max-steps' : [1],
            '-solver' : ['fdlu', 'superlu'],
            '-ts' : ["bdf2", 'midpoint-bdf'],
            '-mesh' : ["ut_square"],
            '-ref' : [1, 2, 3, 4],
            '-decoupled-ms' : [True],
            }

    elif parameter_set == "hlib-speed":
         args_dict = {
            '-binary' : ["./driver/driver"],
            '-driver' : ['llg'],
            '-ts' : ['bdf2'],
            '-mesh' : ["ut_sphere"],
            '-ref' : [1, 2, 3, 4, 5],
            '-tol' : [1e-4],
            '-tmax' : [0.01],
            '-solver' : ['som-gmres'],
            '-prec' : ['som-main-blockut'],
            '-hierarchical-bem' : ['0', '1'],
            }

    else:
        raise NotImplementedError("no parameter set " + str(parameter_set))

    output_root = pjoin('../experiments/parameter_sweeps',
                        '_'.join(parameter_set.split()))


    if cleanup:
        print("Cleaning out", output_root)
        # recursive_check_filenames_rm_safe(output_root)
        shutil.rmtree(output_root)
        os.mkdir(output_root)

    # Make sure the binaries are up to date (if they aren't just the
    # default)
    binaries = args_dict.get('-binary')
    if binaries is not None and not no_build:
        driver_folders = [os.path.abspath(os.path.dirname(d)) for d in binaries]
        for f in driver_folders:
            build_driver(f)

    print("Running parameter sweep with parameter set", parameter_set)
    print("Output is going into", output_root)

    # Run the parameter sweep!
    fun = par(execute_oomph_driver, output_root=output_root)
    mm.parallel_parameter_sweep(fun, args_dict, serial_mode)

    return 0


def main():
    """
    """


    # Parse inputs
    # ============================================================

    parser = argparse.ArgumentParser(description=main.__doc__,

    # Don't mess up my formating in the help message
    formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--debug-mode', action='store_true',
                        help = 'Enable debugging mode (run in serial).')

    parser.add_argument('--parameters', '-p', dest='parameters',
                        help = 'Do a standard parameter sweep with the specified parameter set.')

    parser.add_argument('--clean', action='store_true',
                        help='clean up old results from the target folder')

    parser.add_argument('--no-build', action='store_true',
                        help="Don't rebuild anything")

    args = parser.parse_args()


    if not args.no_build:
        # Make sure micromag library is up to date
        driver_folder = os.path.dirname(mm._DRIVER_PATH)
        library_folder = pjoin(driver_folder, "../../")
        print("Building and installing libraries from", library_folder)
        subp.check_call(['make', '--silent', '--keep-going',
                         'LIBTOOLFLAGS=--silent'], cwd=library_folder)
        subp.check_call(['make', 'install', '--silent', '--keep-going',
                         'LIBTOOLFLAGS=--silent'], cwd=library_folder)

        print("Building driver in", driver_folder)
        subp.check_call(['make', '--silent', '--keep-going',
                         'LIBTOOLFLAGS=--silent'], cwd=driver_folder)


    standard_sweep(args.parameters, args.clean, args.debug_mode,
                   args.no_build)

    return 0


# If this script is run from a shell then run main() and return the result.
if __name__ == "__main__":
    sys.exit(main())
