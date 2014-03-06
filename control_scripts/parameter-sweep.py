#!/usr/bin/env python3

# Future proofing
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

# Imports from main libraries
import sys
import argparse
import os
import os.path
import shutil
import ast

import subprocess as subp
import itertools as it
import functools as ft
import scipy as sp
import multiprocessing as mp


# Imports for specific functions
from functools import partial as par
from os.path import join as pjoin
from glob import glob

# Make sure *this* versions oomphpy is in the path (before any other
# versions in other places)
sys.path.insert(1, pjoin(os.path.dirname(__file__), "../etc"))
import oomphpy
import oomphpy.micromagnetics as mm


def build_driver(folder):
    print("Building in", folder)
    subp.check_call(['make', '--silent', '--keep-going',
                    'LIBTOOLFLAGS=--silent'], cwd=folder)


def main():
    """Run driver multiple times in parallel with arguments specified by a file
    in etc/parameter_sets.
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

    parser.add_argument('--ncores', '-j', '-n',
                        help='Number of processes to run at once')

    parser.add_argument('--no-build', action='store_true',
                        help="Don't rebuild anything")

    args = parser.parse_args()


    if args.ncores is None:
        args.ncores = mp.cpu_count()


    # Find the parameter set
    # ============================================================

    search_root = pjoin(mm.rootdir(), "etc", "parameter_sets")

    # Recurively find files named args.parameters in search_root
    parameter_files = []
    for root, dirs, files in os.walk(search_root, followlinks=True):
        for f in files:
            if f == args.parameters:
                parameter_files.append(pjoin(root, f))
    parameter_file = parameter_files[0]


    # Error check number of files
    if len(parameter_files) > 1:
        sys.stderr.write("Found multiple files named "+ args.parameters + ": "
                         + " ".join(parameter_files) + " .\n")
        return 5
    elif len(parameter_files) == 0:
        sys.stderr.write("Couldn't find a file named "+ args.parameters
                          + " in " + search_root + " .\n")
        return 6


    # Parse parameters file
    # ============================================================

    output_root = pjoin(mm.rootdir(), "experiments", "parameter_sweeps",
                        '_'.join(args.parameters.split()))

    with open(parameter_file, 'r') as pfile:
        args_dict = ast.literal_eval(pfile.read())


    # Make sure we're ready to go
    # ============================================================

    # Maybe build things
    if not args.no_build:

        # Make sure micromag library is up to date
        driver_folder = os.path.dirname(mm.driver_path())
        library_folder = pjoin(driver_folder, "../../")

        print("Building and installing libraries from", library_folder)
        subp.check_call(['make', '--silent', '--keep-going',
                         'LIBTOOLFLAGS=--silent'], cwd=library_folder)
        subp.check_call(['make', 'install', '--silent', '--keep-going',
                         'LIBTOOLFLAGS=--silent'], cwd=library_folder)

        print("Building driver in", driver_folder)
        subp.check_call(['make', '--silent', '--keep-going',
                         'LIBTOOLFLAGS=--silent'], cwd=driver_folder)

        # Make sure the binaries are up to date (if they aren't just the
        # default one).
        binaries = args_dict.get('-binary')
        if binaries is not None:
            driver_folders = [os.path.abspath(os.path.dirname(d)) for d in binaries]
            for f in driver_folders:
                build_driver(f)


    # Remove old stuff if requested
    if args.clean and os.path.isdir(output_root):
        print("Cleaning out", output_root)
        # recursive_check_filenames_rm_safe(output_root)
        shutil.rmtree(output_root)
        os.mkdir(output_root)

    # Copy parameters file to output dir
    os.makedirs(output_root, exist_ok=True)
    shutil.copyfile(parameter_file, pjoin(output_root, "parameter_file"))



    # Run it
    # ============================================================

    print("Running parameter sweep with parameter set", args.parameters)
    print("Output is going into", output_root)
    mm.run_sweep(args_dict, output_root, parallel_sweep=(not args.debug_mode),
                 processes=args.ncores)

    return 0


# If this script is run from a shell then run main() and return the result.
if __name__ == "__main__":
    sys.exit(main())
