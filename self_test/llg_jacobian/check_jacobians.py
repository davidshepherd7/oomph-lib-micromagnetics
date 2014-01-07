#!/usr/bin/env python

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


# Binaries (??ds globals are bad...)
DRIVER = "../../control_scripts/driver/driver"
FPDIFF = "../../../../bin/fpdiff.py"


def dict2argslist(inputdict):
    """For each entry in a dict convert it into a list of the form
       ['-name', 'value', '-name2' 'value2', .... ]
    for use as input to a subprocess command.
    """

    # Support python 2 and 3
    try:
        # python 2 version
        dict_iter = inputdict.iteritems()
    except AttributeError:
        # python 3 version
        dict_iter = inputdict.items()

    processed_kwargs = []
    for key, value in dict_iter:
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

    argslist = dict2argslist(argsdict)

    outdir = argsdict['outdir']
    recreate_dir(outdir)

    validatadir = pjoin('validata', os.path.relpath(outdir, 'Validation'))

    # Run the driver with the given arguments and outputting the Jacobian
    l = [DRIVER, 'llg',
          '-disable-ms',
         '-output-jac', 'always',
         '-solver', 'gmres',
         '-prec', 'amg',
         '-doc-interval', '0',
         ] + argslist
    print("Running", ' '.join(l))
    try:
        subp.check_call(l, stdout=open(pjoin(outdir, 'stdout'), 'w'),
                        stderr=subp.STDOUT)
    except subp.CalledProcessError:
        print("FAILED", ' '.join(l))
        return False

    return True


def check_jacobians(argsdict):

    outdir = argsdict['outdir']
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

    parser.add_argument('--serial', action='store_true')
    parser.add_argument('--fast', action='store_true')
    parser.add_argument('--update-data', action='store_true')

    args = parser.parse_args()


    # Main function
    # ============================================================

    # Build
    print("Building driver")
    subp.check_call(['make'],
                    stdout=open(os.devnull, 'w'),
                    cwd="../../control_scripts/driver/")

    # Set of parameters to test the Jacobians for. Use varying initial m to
    # get mostly non-zeros in J. Only first entry is used for "fast mode".
    jacobian_params = [{'mesh': 'ut_square',
                        'ref': 3,
                        'dt': 1e-4,
                        'tmax': 3e-4,
                        'initm': 'smoothly_varying_5',
                        'outdir': 'Validation/J1'},

                        {'mesh': 'sq_cubeoid',
                        'ref': 2,
                        'dt': 1e-4,
                        'tmax': 2e-4,
                        'initm': 'smoothly_varying_500',
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

        # Make the new Jacobians etc.
        list(Pool().map(generate_jacobians, jacobian_params))

        #  and move them to validata
        files_needed = (glob.glob('Validation/*/jacobian*_1')
                        + glob.glob('Validation/*/residual*_1')
                        + glob.glob('Validation/*/info')
                        + glob.glob('Validation/*/trace')
                        + glob.glob('Validation/*/stdout'))

        for f in files_needed:
            _, d = os.path.split(os.path.dirname(f))
            os.renames(f, pjoin("validata", d, os.path.basename(f)))

        # Zip up all the files
        subp.check_call('gzip validata/*/*', shell=True)

        return 0


    # If its a fast check just do the first one
    if args.fast:
        jacobian_params = [jacobian_params[0]]

    # Otherwise calculate and compare all the Jacobians
    if args.serial:
        results = list(map(jacobian_dump_check, jacobian_params))
    else:
        results = list(Pool().map(jacobian_dump_check, jacobian_params, chunksize=1))

    if all(results):
        return 0
    else:
        return 2


# If this script is run from a shell then run main() and return the result.
if __name__ == "__main__":
    sys.exit(main())
