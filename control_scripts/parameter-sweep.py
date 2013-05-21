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

import itertools as it
import functools as ft
import scipy as sp
import matplotlib.pyplot as plt

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

    # Convert any keyword args into correct format for command line input.
    processed_kwargs = []
    for key, value in args_dict.iteritems():
        if key not in ['mpi_ncores', 'driver', 'outdir']:
            processed_kwargs.append('-'+str(key))
            processed_kwargs.append(str(value))

    # Construct an output directory name based on inputs if not specified.
    if args_dict.get('outdir') is None:
        # Create a has of the inputs and use it to label the folder
        h = hashlib.sha224( ''.join([str(v) for _, v in args_dict.iteritems()]) )
        outdir = ("results_" + h.hexdigest())
    final_outdir = pjoin(output_root, outdir)

    # Make sure the directory is empty and exists
    try:
        os.makedirs(final_outdir)
    except OSError:
        for result_file in glob(pjoin(final_outdir, "*")):
            os.remove(result_file)

    # Construct argument list
    arglist = (['mpirun', '-np', str(args_dict.get('mpi_ncores', 1)),
                str(args_dict['driver']),
                '-outdir', final_outdir]
               + processed_kwargs)

    # Run with specified args in the driver directory, put output (stdout
    # and stderr) into a file.
    with open(pjoin(final_outdir, "stdout"), 'w') as stdout_file:
        err_code = subp.call(arglist,
                             stdout = stdout_file,
                             stderr = subp.STDOUT)

    # Report result
    command = ' '.join(arglist)
    if err_code == 0:
        print(greenColour, command, endColour)
    else:
        # Print failure message to screen
        print('\n', redColour, command , endColour)
        print(redColour, "FAILED with exit code", err_code, "see",
              pjoin(final_outdir, "stdout"), endColour)

        # Print failure message into file in ouput directory
        print('This run failed!', file=open(pjoin(final_outdir, "FAILED"), 'w'))

    return 0


def _apply_to_list(function, list_of_args):
    """Does what it says. Should really be a lambda function but
    multiprocessing requires named functions
    """
    return function(*list_of_args)


def parallel_parameter_sweep(function, parameter_dictionary, serial_mode=False):
    """Run function with all combinations of parameters in parallel using
    all available cores.

    parameter_lists should be a list of lists of parameters,
    """

    import multiprocessing

    # Generate a complete set of combinations of parameters
    parameter_sets = [dict(it.izip(parameter_dictionary, x))
                      for x in it.product(*parameter_dictionary.itervalues())]

    # multiprocessing doesn't include a "starmap", requires all functions
    # to take a single argument. Use a function wrapper to fix this. Also
    # print the list of args while we're in there.

    # For debugging we often need to run in serial (to get useful stack
    # traces).
    if serial_mode:
        results = map(function, parameter_sets)

    else:
        # Run in all parameter sets in parallel
        results = multiprocessing.Pool().map(function, parameter_sets, 1)

    return results


def milan_jacobians(parameter_set, serial_mode=False):

    # Presumably we will always be using the implicit driver...


    # Construct lists of args
    if parameter_set == "initial":
        args_dict = {
            'driver' : "./llg_driver/llg_driver",
            'dt' : [0.1, 0.05, 0.01, 0.001],
            'tmax' : [0.001],
            'tol' : [0.0],
            'refine' : [1, 2, 3, 4, 5],
            'ts' : ["bdf2"],
            'initm' : ['smoothly_varying'],
            'happ' : ['x', 'y', 'z'],
            'mesh' : ['sq_square', 'ut_square'],
            'output_jac' : ['at_end']
            }
        # more parameter sets go here
    else:
        raise NotImplementedError

    # Fill in some function args that should always be the same.
    fun = par(execute_oomph_driver,
              output_root='../experiments/jacobian_sweeps')

    # Run the parameter sweep!
    parallel_parameter_sweep(fun, args_dict, serial_mode)

    return 0

def standard_sweep(parameter_set, serial_mode=False):

    # Construct lists of args
    if parameter_set == '0':
        args_dict = {
            'driver' : ["./llg_driver/llg_driver"],
            'dt' : [0.1, 0.05, 0.01, 0.001],
            'tmax' : [6.0],
            'tol' : [0.0],
            'refine' : [1, 2, 3, 4, 5],
            'ts' : ["bdf2", 'midpoint'],
            'initm' : ['z', 'smoothly_varying'],
            'happ' : ['minus_z'],
            'mesh' : ['sq_square', 'ut_square'],
            'solver' : ['superlu'],
            }

    elif parameter_set == '1':
        args_dict = {
            'driver' : ["./llg_driver/llg_driver"],
            'dt' : [0.1, 0.01],
            'tmax' : [2.0],
            'tol' : [0.0],
            'refine' : [1, 5],
            'ts' : ['midpoint'],
            'initm' : ['z', 'smoothly_varying'],
            'happ' : ['minus_z'],
            'mesh' : ['ut_square'],
            'solver' : ['superlu'],
            }

    elif parameter_set == '2':
        args_dict = {
            'driver' : ["./llg_driver/llg_driver"],
            'dt' : [1e-6],
            'tmax' : [2.0],
            'tol' : [1e-3, 1e-4, 1e-5],
            'refine' : [1, 2, 3],
            'ts' : ['bdf2', 'midpoint'],
            'initm' : ['z', 'smoothly_varying'],
            'happ' : ['minus_z'],
            'mesh' : ['ut_square'],
            'solver' : ['superlu'],
            }

    elif parameter_set == '3':
        args_dict = {
            'driver' : ["./llg_driver/llg_driver"],
            'dt' : [1e-6],
            'tmax' : [2.0],
            'tol' : [1e-3, 1e-4, 1e-5],
            'refine' : [1, 2, 3],
            'ts' : ['bdf2', 'midpoint'],
            'initm' : ['z'],
            'happ' : ['minus_z'],
            'mesh' : ['ut_square', 'sq_square'],
            'solver' : ['superlu'],
            }

    elif parameter_set == '4':
        args_dict = {
            'driver' : ["./llg_driver/llg_driver"],
            'dt' : [1e-6],
            'tmax' : [2.0],
            'tol' : [1e-3, 1e-4, 1e-5],
            'refine' : [1, 2, 3],
            'ts' : ['bdf2', 'midpoint'],
            'initm' : ['z'],
            'happ' : ['minus_z'],
            'mesh' : ['ut_sphere'],
            'solver' : ['superlu'],
            }

    elif parameter_set == 'cubeoid':
        args_dict = {
            'driver' : ["./semi_implicit_mm_driver/semi_implicit_mm_driver"],
            'dt' : [1e-6],
            'tmax' : [2.0],
            'tol' : [1e-3, 1e-4, 1e-5],
            'refine' : [1, 2, 3],
            'ts' : ['bdf2', 'midpoint'],
            'initm' : ['z'],
            'happ' : ['minus_z'],
            'mesh' : ['ut_cubeoid', 'st_cubeoid'],
            'solver' : ['superlu'],
            }

    elif parameter_set == 'nmag_cubeoid':
        args_dict = {
            'driver' : ["./semi_implicit_mm_driver/semi_implicit_mm_driver"],
            'dt' : [1e-4],
            'tmax' : [60.0],
            'tol' : [1e-3, 1e-4, 1e-5],
            'refine' : [4],
            'ts' : ['midpoint', 'bdf2'],
            'initm' : ['xz'],
            'happ' : ['zero'],
            'mesh' : ['sq_cubeoid'],
            'solver' : ['gmres-amg'],
            }

    elif parameter_set == 'const_dt_nmag_cubeoid':
        args_dict = {
            'driver' : ["./semi_implicit_mm_driver/semi_implicit_mm_driver"],
            'dt' : [1e-4, 5e-5],
            'tmax' : [30.0],
            'tol' : [0.0],
            'refine' : [4],
            'ts' : ['midpoint', 'bdf2'],
            'initm' : ['xz'],
            'happ' : ['zero'],
            'mesh' : ['sq_cubeoid'],
            'solver' : ['gmres-amg'],
            }

    elif parameter_set == 'unsteady_heat_midpoint_vs_bdf2':
        args_dict = {
            'driver' : ["./unsteady_heat_driver/unsteady_heat_driver"],
            'dt' : [1e-4],
            'tmax' : [10.0],
            'tol' : [1e-3, 1e-4, 1e-5],
            'ts' : ['midpoint', 'bdf2'],
            }

    elif parameter_set == 'check_semi_impl':
        args_dict = {
            'driver' : ["./semi_implicit_mm_driver/semi_implicit_mm_driver"],
            'dt' : [1e-6],
            'tmax' : [2.0],
            'tol' : [1e-3],
            'refine' : [1,2],
            'ts' : ['midpoint'],
            'initm' : ['z', 'smoothly_varying'],
            'happ' : ['minus_z'],
            'mesh' : ['ut_sphere', 'sq_square'],
            'solver' : ['superlu'],
            }

    elif parameter_set == 'cubeoid-timestep-newton-convergence':
        args_dict = {
            'driver' : ["./llg_driver/llg_driver"],
            'dt' : [1e-4, 1e-5, 1e-6, 1e-7, 1e-8],
            'tmax' : [1e-8],
            'tol' : [0.0],
            'refine' : [3],
            'ts' : ['bdf2'],
            'initm' : ['z'],
            'happ' : ['minus_z'],
            'mesh' : ['ut_cubeoid', 'sq_cubeoid', 'st_cubeoid'],
            'solver' : ['superlu'],
            }

    else:
        raise NotImplementedError("no parameter set " + str(parameter_set))

    output_root = pjoin('../experiments/parameter_sweeps',
                        '_'.join(parameter_set.split()))

    print("Running parameter sweep with parameter set", parameter_set)
    print("Output is going into", output_root)

    # Run the parameter sweep!
    fun = par(execute_oomph_driver, output_root=output_root)
    parallel_parameter_sweep(fun, args_dict, serial_mode)

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

    parser.add_argument('--parameters', '-p', dest='parameters',
                        help = 'Do a standard parameter sweep with the specified parameter set.')

    args = parser.parse_args()

    # Main function
    # ============================================================

    # Make sure the driver binary is up to date
    print("Building in ./llg_driver folder.")
    subp.check_call(['make', '--silent', '--keep-going',
                     'LIBTOOLFLAGS=--silent'], cwd = "./llg_driver")

    # print("Building in ./semi_implicit_mm_driver folder.")
    # subp.check_call(['make', '--silent', '--keep-going',
    #                  'LIBTOOLFLAGS=--silent'], cwd = "./semi_implicit_mm_driver")

    print("Building in ./unsteady_heat_driver folder.")
    subp.check_call(['make', '--silent', '--keep-going',
                     'LIBTOOLFLAGS=--silent'], cwd = "./unsteady_heat_driver")

    if args.parameters is not None:
        standard_sweep(args.parameters, args.debug_mode)

    elif args.j_parameter_set is not None:
        print("Running Jacobian generation parameter sweep",
              "with parameter set", args.j_parameter_set)
        milan_jacobians(args.j_parameter_set,args.debug_mode)

    return 0


# If this script is run from a shell then run main() and return the result.
if __name__ == "__main__":
    sys.exit(main())
