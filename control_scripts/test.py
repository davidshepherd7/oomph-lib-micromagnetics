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

import itertools as it
import functools as ft
import scipy as sp
import matplotlib.pyplot as plt

# Imports for specific functions
from functools import partial as par
from os.path import join as pjoin
from glob import glob


# TODO:

# On error put stdout file somewhere useful.

def execute_oomph_driver(mpi_ncores, driver, dt, tmax, tol, refinement,
                         outdir, timestepper, initial_m, applied_field, mesh,
                         output_root = "./",
                         **kwargs):

    # Convert keyword args into correct format for command line input.
    processed_kwargs = []
    for key, value in kwargs.iteritems():
        processed_kwargs.append('-'+str(key))
        processed_kwargs.append(str(value))

    # Construct an output directory name based on inputs if not specified.
    if outdir is None:
        outdir = ("results_" + str(dt) + "_" + str(tol) + "_" + str(refinement)
                  + "_" + timestepper + "_" + applied_field + "_" + mesh
                  + "_" + initial_m)
    final_outdir = os.path.join(output_root, outdir, "")


    # Make sure the directory is empty and exists
    try:
        os.makedirs(final_outdir)
    except OSError:
        for result_file in glob(pjoin(outdir,"*")):
            os.remove(result_file)

    # Run with specified args, put output (stdout and stderr) into a file.
    with open(pjoin(final_outdir, "stdout"), 'w') as stdout_file:
        arglist = ['mpirun', '-np', str(mpi_ncores),
                   driver,
                   '-dt', str(dt),
                   "-tmax", str(tmax),
                   "-tol", str(tol),
                   "-ref", str(refinement),
                   "-outdir", final_outdir,
                   "-ts", timestepper,
                   "-initm", initial_m,
                   "-happ", applied_field,
                   "-mesh", mesh] + processed_kwargs

        err_code = subp.call(arglist, stdout = stdout_file,
                             stderr = subp.STDOUT)

    if err_code != 0:
        print(final_outdir, "FAILED with exit code", err_code)

    return 0


def _apply_to_list_and_print_args(function, list_of_args):
    """Does what it says. Should really be a lambda function but
    multiprocessing requires named functions
    """
    print(list_of_args)
    return function(*list_of_args)


def parallel_parameter_sweep(function, parameter_lists, serial_mode=False):
    """Run function with all combinations of parameters in parallel using
    all available cores.

    parameter_lists should be a list of lists of parameters,
    """

    import multiprocessing

    # Generate a complete set of combinations of parameters
    parameter_sets = it.product(*parameter_lists)

    # multiprocessing doesn't include a "starmap", requires all functions
    # to take a single argument. Use a function wrapper to fix this. Also
    # print the list of args while we're in there.
    wrapped_function = par(_apply_to_list_and_print_args, function)

    # For debugging we often need to run in serial (to get useful stack
    # traces).
    if serial_mode:
        results_iterator = it.imap(wrapped_function, parameter_sets)

    else:
        # Run in all parameter sets in parallel
        pool = multiprocessing.Pool()
        results_iterator = pool.imap_unordered(wrapped_function, parameter_sets)
        pool.close()

        # wait for everything to finish
        pool.join()

    # Force immediate evaluation of the functions (imap is lazily
    # evaluated)
    results = list(results_iterator)

    return results


def milan_jacobians(parameter_set, serial_mode=False):

    # Construct lists of args
    if parameter_set == "initial":
        dts = [0.1, 0.05, 0.01, 0.001]
        tmaxs = [0.001]
        tols = [0.0]
        refines = [1, 2, 3, 4, 5]
        outdirs = [None]
        timesteppers = ["bdf2"]
        initial_ms = ['smoothly_varying']
        fields = ['x', 'y', 'z']
        meshes = ['sq_square', 'ut_square']

        # more parameter sets go here
    else:
        raise NotImplementedError

    arg_list = [dts, tmaxs, tols, refines, outdirs, timesteppers,
                initial_ms, fields, meshes]

    # Fill in some function args that should always be the same.
    fun = par(execute_oomph_driver, 1, "./driver/driver",
              output_root='../experiments/jacobian_sweeps',
              output_jac='at_end')

    # Run the parameter sweep!
    parallel_parameter_sweep(fun, arg_list, serial_mode)

    return 0

def midpoint_comparisons(parameter_set, serial_mode=False):

    # Construct lists of args
    dts = [0.1, 0.05, 0.01, 0.001]
    tmaxs = [6.0]
    tols = [0.0]
    refines = [1, 2, 3, 4, 5]
    outdirs = [None]
    timesteppers = ["bdf2", 'midpoint']
    initial_ms = ['z', 'smoothly_varying']
    fields = ['minus_z']
    meshes = ['sq_square', 'ut_square']

    arg_list = [dts, tmaxs, tols, refines, outdirs, timesteppers,
                initial_ms, fields, meshes]

    # Fill in some function args that should always be the same.
    fun = par(execute_oomph_driver, 1, "./driver/driver",
          output_root='../experiments/midpoint_sweeps')

        # Run the parameter sweep!
    parallel_parameter_sweep(fun, arg_list, serial_mode)

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

    parser.add_argument('--jacobians', dest='j_parameter_set',
                        help = 'Do a parameter sweep just dumping jacobians.')
    parser.add_argument('--midpoint', dest='midpoint_parameter_set',
                        help = 'Do a parameter sweep comparing midpoint and BDF2')

    args = parser.parse_args()

    # Main function
    # ============================================================

    # Make sure the driver binary is up to date
    print("Building in ./driver folder.")
    subp.check_call(['make', '--silent', '--keep-going',
                     'LIBTOOLFLAGS=--silent'], cwd = "./driver")


    if args.j_parameter_set is not None:
        milan_jacobians(args.j_parameter_set, args.debug_mode)

    if args.midpoint_parameter_set is not None:
        midpoint_comparisons(args.midpoint_parameter_set, args.debug_mode)

    return 0


# If this script is run from a shell then run main() and return the result.
if __name__ == "__main__":
    sys.exit(main())
