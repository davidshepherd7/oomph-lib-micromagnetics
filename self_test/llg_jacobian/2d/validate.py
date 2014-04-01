#!/usr/bin/env python3

# Future proofing
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals


# Imports from main libraries
import subprocess as subp
from multiprocessing import Pool
import multiprocessing
import sys
import argparse
import os
import os.path
import glob
import random

import itertools as it
import functools as ft
import scipy as sp

# Imports for specific functions
from functools import partial as pt
from os.path import join as pjoin

# Make sure *this* versions oomphpy is in the path (before any other
# versions in other places)
sys.path.insert(1, pjoin(os.path.dirname(__file__), "../../../etc"))
import oomphpy
import oomphpy.micromagnetics as mm


# Binaries (??ds globals are bad...)
FPDIFF = "../../../../../bin/fpdiff.py"


def generate_jacobians(argdict):

    # Add fixed values to argument dictionary
    dict2 = {'-driver' : 'llg',
            '-ms-method' : 'disabled',
            '-output-jac' : 'always',
            '-solver' : 'gmres',
            '-prec' : 'amg',
            '-doc-interval' : '0'}
    argdict.update(dict2)

    # Get the arguments as a list
    arglist, binary_path, _ = mm.argdict2list(argdict)

    # Make sure the output dir exists and is clean
    outdir = argdict['-outdir']
    mm.cleandir(outdir)

    # Run the driver with the given arguments
    err_code = mm.run_driver(arglist, outdir)

    # Return sucess/failure
    if err_code != 0:
        print("FAILED", ' '.join(arglist))
        return False
    else:
        return True


def check_jacobians(argsdict):

    outdir = argsdict['-outdir']
    validatadir = pjoin('validata', os.path.relpath(outdir, 'Validation'))

    # Compare all Jacobians using fpdiff
    compare_files = (glob.glob(pjoin(validatadir, 'jacobian*'))
                     + glob.glob(pjoin(validatadir, 'residual*')))

    success = True
    with open(pjoin(outdir, 'check_trace'), 'a') as tracefile:
        for stored_file in compare_files:

            # Calculate filenames: strip .gz and change dir to Validation
            basename = os.path.basename(os.path.splitext(stored_file)[0])
            computed_file = pjoin(outdir,basename)

            try:
                print("Comparing", stored_file, computed_file, end=' ')
                subp.check_call([FPDIFF, stored_file, computed_file],
                                stdout=tracefile,
                                stderr=tracefile)
            except subp.CalledProcessError as e:
                if e.returncode == 2:
                    print("FAILED")
                    success = False
                else:
                    raise
            else:
                print() # Finsh the line
    return success


def jacobian_dump_check(argsdict):
    return generate_jacobians(argsdict) and check_jacobians(argsdict)


def main():
    """

    """


    # Parse inputs
    # ============================================================

    parser = argparse.ArgumentParser(description=main.__doc__,

    # Don't mess up my formating in the help message
    formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--update-data', action='store_true')

    args = parser.parse_args()


    # Main function
    # ============================================================

    # # Build
    # print("Building driver")
    # subp.check_call(['make'],
    #                 stdout=open(os.devnull, 'w'),
    #                 cwd="../../control_scripts/driver/")

    # Set of parameters to test the Jacobians for. Use varying initial m to
    # get mostly non-zeros in J.
    jacobian_params = {'-mesh': 'sq_square',
                        '-ref': 3,
                        '-dt': 1e-4,
                        '-tmax': 3e-4,
                        '-initial-m': 'smoothly_varying_5',
                        '-outdir': 'Validation'}


    # If we are updating the data just do that then exit
    if args.update_data:

        # Move old Jacobians
        old_jac_folder = 'validata.bak.' + str(random.random())
        #??ds use something based on git hash instead
        print("Replacing Jacobian data with newly generated Jacobians!")
        try:
            os.renames('validata', old_jac_folder)
        except OSError:
            print("No old Jacobians to move.")
        else:
            print("Old Jacobians are moving to", old_jac_folder)

        # Make the new Jacobians etc.
        generate_jacobians(jacobian_params)

        #  and move them to validata
        files_needed = (glob.glob('Validation/jacobian*_1')
                        + glob.glob('Validation/residual*_1')
                        + ['Validation/info', 'Validation/trace',
                            'Validation/stdout'])

        for f in files_needed:
            os.renames(f, pjoin("validata", os.path.basename(f)))

        # Zip up all the files
        subp.check_call('gzip validata/*', shell=True)

        return 0


    result = jacobian_dump_check(jacobian_params)

    if result:
        return 0
    else:
        return 2


# If this script is run from a shell then run main() and return the result.
if __name__ == "__main__":
    sys.exit(main())
