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


# possible improvements:

# On error put stdout file somewhere useful?


greenColour = '\033[01;32m'
redColour = '\033[01;31m'
endColour = '\033[0m'

def execute_oomph_driver(mpi_ncores, driver, dt, tmax, tol, refinement,
                         outdir, timestepper, initial_m, applied_field, mesh,
                         solver,
                         output_root = "./",
                         **kwargs):

    # Convert any keyword args into correct format for command line input.
    processed_kwargs = []
    for key, value in kwargs.iteritems():
        processed_kwargs.append('-'+str(key))
        processed_kwargs.append(str(value))

    # Construct an output directory name based on inputs if not specified.
    if outdir is None:
        outdir = ("results_" + str(dt) + "_" + str(tol) + "_" + str(refinement)
                  + "_" + timestepper + "_" + applied_field + "_" + mesh
                  + "_" + initial_m)
    final_outdir = pjoin(output_root, outdir)

    # Make sure the directory is empty and exists
    try:
        os.makedirs(final_outdir)
    except OSError:
        for result_file in glob(pjoin(final_outdir, "*")):
            os.remove(result_file)

    # Run with specified args in the driver directory, put output (stdout
    # and stderr) into a file.
    with open(pjoin(final_outdir, "stdout"), 'w') as stdout_file:
        arglist = (['mpirun', '-np', str(mpi_ncores),
                   driver,
                   '-dt', str(dt),
                   "-tmax", str(tmax),
                   "-tol", str(tol),
                   "-ref", str(refinement),
                   "-outdir", final_outdir,
                   "-ts", timestepper,
                   "-initm", initial_m,
                   "-happ", applied_field,
                   "-mesh", mesh,
                   "-solver", solver]
                   + processed_kwargs)

        err_code = subp.call(arglist,
                             stdout = stdout_file,
                             stderr = subp.STDOUT)

    if err_code == 0:
        print(greenColour, [os.path.basename(driver), dt, tmax, tol, refinement,
                            timestepper, initial_m, applied_field, mesh],
                            endColour)
    else:
        # Print failure message to screen
        print('\n', redColour, ' '.join(arglist), endColour)
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
    wrapped_function = par(_apply_to_list, function)

    # For debugging we often need to run in serial (to get useful stack
    # traces).
    if serial_mode:
        results = map(wrapped_function, parameter_sets)

    else:
        # Run in all parameter sets in parallel
        pool = multiprocessing.Pool()
        results = pool.map(wrapped_function, parameter_sets, 1)
        pool.close()

    return results


def milan_jacobians(parameter_set, serial_mode=False):

    # Presumably we will always be using the implicit driver...
    rel_driver_path = "./llg_driver/llg_driver"

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
    fun = par(execute_oomph_driver, 1,
              rel_driver_path,
              output_root='../experiments/jacobian_sweeps',
              output_jac='at_end')

    # Run the parameter sweep!
    parallel_parameter_sweep(fun, arg_list, serial_mode)

    return 0

def standard_sweep(parameter_set, serial_mode=False):

    # Construct lists of args
    if parameter_set == '0':
        rel_driver_paths = ["./llg_driver/llg_driver"]
        dts = [0.1, 0.05, 0.01, 0.001]
        tmaxs = [6.0]
        tols = [0.0]
        refines = [1, 2, 3, 4, 5]
        outdirs = [None]
        timesteppers = ["bdf2", 'midpoint']
        initial_ms = ['z', 'smoothly_varying']
        fields = ['minus_z']
        meshes = ['sq_square', 'ut_square']
        solvers = ['superlu']

    elif parameter_set == '1':
        rel_driver_paths = ["./llg_driver/llg_driver"]
        dts = [0.1, 0.01]
        tmaxs = [2.0]
        tols = [0.0]
        refines = [1, 5]
        outdirs = [None]
        timesteppers = ['midpoint']
        initial_ms = ['z', 'smoothly_varying']
        fields = ['minus_z']
        meshes = ['ut_square']
        solvers = ['superlu']

    elif parameter_set == '2':
        rel_driver_paths = ["./llg_driver/llg_driver"]
        dts = [1e-6]
        tmaxs = [2.0]
        tols = [1e-3, 1e-4, 1e-5]
        refines = [1, 2, 3]
        outdirs = [None]
        timesteppers = ['bdf2', 'midpoint']
        initial_ms = ['z', 'smoothly_varying']
        fields = ['minus_z']
        meshes = ['ut_square']
        solvers = ['superlu']

    elif parameter_set == '3':
        rel_driver_paths = ["./llg_driver/llg_driver"]
        dts = [1e-6]
        tmaxs = [2.0]
        tols = [1e-3, 1e-4, 1e-5]
        refines = [1, 2, 3]
        outdirs = [None]
        timesteppers = ['bdf2', 'midpoint']
        initial_ms = ['z']
        fields = ['minus_z']
        meshes = ['ut_square', 'sq_square']
        solvers = ['superlu']

    elif parameter_set == '4':
        rel_driver_paths = ["./llg_driver/llg_driver"]
        dts = [1e-6]
        tmaxs = [2.0]
        tols = [1e-3, 1e-4, 1e-5]
        refines = [1, 2, 3]
        outdirs = [None]
        timesteppers = ['bdf2', 'midpoint']
        initial_ms = ['z']
        fields = ['minus_z']
        meshes = ['ut_sphere']
        solvers = ['superlu']

    elif parameter_set == 'cubeoid':
        rel_driver_paths = ["./semi_implicit_mm_driver/semi_implicit_mm_driver"]
        dts = [1e-6]
        tmaxs = [2.0]
        tols = [1e-3, 1e-4, 1e-5]
        refines = [1, 2, 3]
        outdirs = [None]
        timesteppers = ['bdf2', 'midpoint']
        initial_ms = ['z']
        fields = ['minus_z']
        meshes = ['ut_cubeoid', 'st_cubeoid']
        solvers = ['superlu']

    elif parameter_set == 'nmag_cubeoid':
        rel_driver_paths = ["./semi_implicit_mm_driver/semi_implicit_mm_driver"]
        dts = [1e-6]
        tmaxs = [5.0]
        tols = [1e-3]
        refines = [4]
        outdirs = [None]
        timesteppers = ['midpoint']
        initial_ms = ['xz']
        fields = ['zero']
        meshes = ['sq_cubeoid']
        solvers = ['gmres-amg']

        # ./semi_implicit_mm_driver/semi_implicit_mm_driver -dt 1e-06 -tmax 5.0 -tol 0.001 -ref 5 -outdir ../experiments/parameter_sweeps/nmag_cubeoid_final -initm xz -happ zero -mesh sq_cubeoid -solver gmres-amg

    elif parameter_set == 'check_semi_impl':
        rel_driver_paths = ["./semi_implicit_mm_driver/semi_implicit_mm_driver"]
        dts = [1e-6]
        tmaxs = [2.0]
        tols = [1e-3]
        refines = [1,2]
        outdirs = [None]
        timesteppers = ['midpoint']
        initial_ms = ['z', 'smoothly_varying']
        fields = ['minus_z']
        meshes = ['ut_sphere', 'sq_square']
        solvers = ['superlu']

    elif parameter_set == 'cubeoid-timestep-newton-convergence':
        rel_driver_paths = ["./llg_driver/llg_driver"]
        dts = [1e-4, 1e-5, 1e-6, 1e-7, 1e-8]
        tmaxs = [1e-8]
        tols = [0.0]
        refines = [3]
        outdirs = [None]
        timesteppers = ['bdf2']
        initial_ms = ['z']
        fields = ['minus_z']
        meshes = ['ut_cubeoid', 'sq_cubeoid', 'st_cubeoid']
        solvers = ['superlu']

    else:
        raise NotImplementedError("no parameter set " + str(parameter_set))

    output_root = pjoin('../experiments/parameter_sweeps',
                        '_'.join(parameter_set.split()))

    print("Running parameter sweep with parameter set", parameter_set)
    print("Output is going into", output_root)

    # Run the parameter sweep!
    fun = par(execute_oomph_driver, 1, output_root=output_root)
    arg_list = [rel_driver_paths, dts, tmaxs, tols, refines, outdirs,
                timesteppers, initial_ms, fields, meshes, solvers]
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

    parser.add_argument('--parameters', '-p', dest='parameters',
                        help = 'Do a standard parameter sweep with the specified parameter set.')

    args = parser.parse_args()

    # Main function
    # ============================================================

    # Make sure the driver binary is up to date
    print("Building in ./llg_driver folder.")
    subp.check_call(['make', '--silent', '--keep-going',
                     'LIBTOOLFLAGS=--silent'], cwd = "./llg_driver")

    print("Building in ./semi_implicit_mm_driver folder.")
    subp.check_call(['make', '--silent', '--keep-going',
                     'LIBTOOLFLAGS=--silent'], cwd = "./semi_implicit_mm_driver")

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
