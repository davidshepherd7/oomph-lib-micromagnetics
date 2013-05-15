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
import glob
import random

import itertools as it
import functools as ft
import scipy as sp
import matplotlib.pyplot as plt

# Imports for specific functions
from functools import partial as pt
from os.path import join as pjoin


# Binaries (??ds globals are bad...)
DRIVER = "../../control_scripts/llg_driver/llg_driver"
FPDIFF = "../../../../bin/fpdiff.py"


def dict2argslist(inputdict):
    """For each entry in a dict convert it into a list of the form
       ['-name', 'value', '-name2' 'value2', .... ]
    for use as input to a subprocess command.
    """

    processed_kwargs = []
    for key, value in inputdict.iteritems():
        processed_kwargs.append('-'+str(key))
        processed_kwargs.append(str(value))

    return processed_kwargs


def recreate_dir(dirname):
    """Remove dir and remake it. Ignore errors from dir not existing before
    removal. Automatically create parents if needed.
    """
    import shutil
    shutil.rmtree(dirname, ignore_errors=True)
    os.makedirs(dirname)


def generate_jacobians(argsdict):

    print("Running", argsdict)

    argslist = dict2argslist(argsdict)

    outdir = argsdict['outdir']
    recreate_dir(outdir)

    validatadir = pjoin('validata', os.path.relpath(outdir, 'Validation'))

    # Run the driver with the given arguments and outputting the Jacobian
    subp.check_call([DRIVER, '-output_jac', 'always',
                     '-solver', 'gmres-amg'] + argslist,
                     stdout=open(pjoin(outdir, 'stdout'), 'w')
        )

    # Remove time data from the info file
    computed_info_file = pjoin(outdir, 'info')
    data="".join(open(computed_info_file).readlines()[5:-1])
    open(computed_info_file, "wb").write(data)

    return


def check_jacobians(argsdict):

    outdir = argsdict['outdir']
    validatadir = pjoin('validata', os.path.relpath(outdir, 'Validation'))

    # Compare the info file and all Jacobians using fpdiff
    compare_files = ([pjoin(validatadir, 'info.gz')]
                     + glob.glob(pjoin(validatadir, 'jacobian*')))

    with open(pjoin(outdir, 'check_trace'), 'a') as tracefile:
        for stored_file in compare_files:

            # Calculate filenames: strip .gz and change dir to Validation
            basename = os.path.basename(os.path.splitext(stored_file)[0])
            computed_file = pjoin(outdir,basename)

            print("Comparing", stored_file, computed_file)
            subp.check_call([FPDIFF, stored_file, computed_file],
                             stdout=tracefile)
    return


def jacobian_dump_check(argsdict):
    generate_jacobians(argsdict)
    check_jacobians(argsdict)
    return


def main():
    """

    """


    # Parse inputs
    # ============================================================

    parser = argparse.ArgumentParser(description=main.__doc__,

    # Don't mess up my formating in the help message
    formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--serial', action='store_true')
    parser.add_argument('--fast', action='store_true')
    parser.add_argument('--update-data', action='store_true')

    args = parser.parse_args()


    # Main function
    # ============================================================

    # Build
    print("Building llg_driver")
    subp.check_call(['make'],
                    stdout=open(os.devnull, 'w'),
                    cwd="../../control_scripts/llg_driver/")

    # Set of parameters to test the Jacobians for. Use varying initial m to
    # get mostly non-zeros in J. Only first entry is used for "fast mode".
    jacobian_params = [{'mesh': 'ut_square',
                        'ref': 3,
                        'dt': 1e-4,
                        'tmax': 3e-4,
                        'initm': 'smoothly_varying',
                        'outdir': 'Validation/J1'},

                        {'mesh': 'sq_cubeoid',
                        'ref': 2,
                        'dt': 1e-4,
                        'tmax': 2e-4,
                        'initm': 'smoothly_varying',
                        'outdir': 'Validation/J2'},

                        ]

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

        # Make the new Jacobians and move them to validata
        list(Pool().map(generate_jacobians, jacobian_params))
        os.renames('Validation', 'validata')

        # Zip up all the files
        subp.check_call('gzip validata/*/*', shell=True)
        return 0

    # If its a fast check just do the first one
    if args.fast:
        jacobian_params = [jacobian_params[0]]

    # Otherwise calculate and compare all the Jacobians
    if args.serial:
        list(map(jacobian_dump_check, jacobian_params))
    else:
        list(Pool().map(jacobian_dump_check, jacobian_params, chunksize=1))


    # If it gets to here then everything worked!
    return 0


# If this script is run from a shell then run main() and return the result.
if __name__ == "__main__":
    sys.exit(main())
