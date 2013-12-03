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

import itertools as it
import functools as ft
import scipy as sp

# Imports for specific functions
from functools import partial as par
from os.path import join as pjoin
from glob import glob


# possible improvements:

# On error put stdout file somewhere useful?


greenColour = '\033[01;32m'
redColour = '\033[01;31m'
endColour = '\033[0m'


# def recursive_check_filenames_rm_safe(root_dir):
#     """Check that all files inside a directory are only output from my drivers,
#     and so safe to be deleted.
#     """
#     whitelist = ["info", "trace", "stdout", "soln.pvd", "run_script.sh",
#                  ".dummy_check.dat", "FAILED",
#                  "info.gz", "trace.gz", "stdout.gz", ".dummy_check.dat.gz",
#                  "FAILED.gz", "averages", "errors", "field_averages",
#                  "solndatfiles.tar.bz2"]

#     regex_whitelist = map(re.compile,
#                           [r"\A(soln|phione|field)(\d*)\.(dat|dat\.gz|vtu)$",
#                            r"jacobian.*", r"residual.*", r".*~"])

#     def safe(f):
#         safe_wl = f in whitelist
#         safe_re = any([(r.match(f) is not None) for r in regex_whitelist])
#         print(safe_wl, safe_re)
#         return safe_wl or safe_re

#     for root, dirs, files in os.walk(root_dir):
#         for f in files:
#             if not safe(f):
#                 print(f)
#                 raise Exception("Tried to delete non whitelisted file: \""+f+"\"")


def execute_oomph_driver(args_dict, output_root):

    # Convert any keyword args into correct format for command line input.
    processed_kwargs = []
    for key, value in args_dict.items():
        if key not in ['mpi_ncores', 'driver', 'outdir', 'implicit-ms']:
            processed_kwargs.append('-'+str(key))
            processed_kwargs.append(str(value))

    # Handle implicit ms: append to args list if true in dict and driver is llg
    # driver.
    implicit_ms = args_dict.get('implicit-ms')
    if (implicit_ms is not None) and args_dict['driver'] == "./llg_driver/llg_driver":
        if implicit_ms:
            processed_kwargs.append('-implicit-ms')


    # Construct an output directory name based on inputs if not specified.
    if args_dict.get('outdir') is None:
        # Create a hash of the inputs and use it to label the folder
        h = hashlib.sha224( ''.join([str(v) for _, v in args_dict.items()]).encode())
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

    # Write the command used to a file
    with open(pjoin(final_outdir, "run_script.sh"), 'w') as script_file:
        script_file.write("#!/bin/sh\n")
        script_file.write("# Command used in run\n")
        script_file.write(' '.join(arglist))
        script_file.write("\n")


    # Run with specified args in the driver directory, put output (stdout
    # and stderr) into a file.
    with open(pjoin(final_outdir, "stdout"), 'w') as stdout_file:
        err_code = subp.call(arglist,
                             stdout = stdout_file,
                             stderr = subp.STDOUT)

    # Compress the output data files (but not info or trace).
    for datfilename in glob(pjoin(final_outdir, "*.dat")):
        subp.check_call(['gzip', datfilename])

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


def build(folder):
    print("Building in", folder)
    subp.check_call(['make', '--silent', '--keep-going',
                    'LIBTOOLFLAGS=--silent'], cwd=folder)


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
    parameter_sets = [dict(zip(parameter_dictionary, x))
                      for x in it.product(*parameter_dictionary.values())]

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
            'ref' : [1, 2, 3, 4, 5],
            'ts' : ["bdf2"],
            'initm' : ['smoothly_varying'],
            'happ' : ['x', 'y', 'z'],
            'mesh' : ['sq_square', 'ut_square'],
            'output-jac' : ['at_end']
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

def standard_sweep(parameter_set, cleanup, serial_mode=False):

    if parameter_set == 'nmag_cubeoid':
        args_dict = {
            'driver' : ["./semi_implicit_mm_driver/semi_implicit_mm_driver"],
            'dt' : [1e-4],
            'tmax' : [60.0],
            'tol' : [1e-3, 1e-4, 1e-5],
            'ref' : [4],
            'ts' : ['midpoint', 'bdf2'],
            'initm' : ['xz'],
            'happ' : ['zero'],
            'mesh' : ['sq_cubeoid'],
            'solver' : ['gmres'],
            'preconditioner' : ['amg']
            }

    elif parameter_set == 'const_dt_nmag_cubeoid':
        args_dict = {
            'driver' : ["./semi_implicit_mm_driver/semi_implicit_mm_driver"],
            'dt' : [1e-2, 5e-3],
            'tmax' : [20.0],
            'tol' : [0.0],
            'ref' : [2, 3, 4],
            'ts' : ['midpoint', 'bdf2'],
            'initm' : ['xz'],
            'happ' : ['zero'],
            'mesh' : ['sq_cubeoid'],
            'solver' : ['gmres'],
            'preconditioner' : ['amg']
            }

    elif parameter_set == 'nmag_cubeoid_llg_prec':
        args_dict = {
            'driver' : ["./semi_implicit_mm_driver/semi_implicit_mm_driver"],
            'dt' : [5e-4],
            'tmax' : [20.0],
            'tol' : [1e-4],
            'ref' : [4],
            'ts' : ['bdf2'],
            'initm' : ['xz'],
            'happ' : ['zero'],
            'mag-params' :["simple-llg"],
            'mesh' : ['sq_cubeoid'],
            'solver' : ['gmres'],
            'preconditioner' : ['blockllg-uppertriangular-blockexact-xy-exact',
                                'blockllg-uppertriangular-blockexact-xz-exact',
                                'blockllg-uppertriangular-blockexact-yz-exact',
                                ]
            }

    elif parameter_set == 'script_test':
        args_dict = {
            'driver' : ["./llg_driver/llg_driver"],
            'dt' : [1e-4],
            'tmax' : [1.0],
            'tol' : [1e-3],
            'ref' : [2],
            }

    elif parameter_set == 'unsteady_heat_midpoint_vs_bdf2':
        args_dict = {
            'driver' : ["./unsteady_heat_driver/unsteady_heat_driver"],
            'dt' : [1e-4],
            'tmax' : [10.0],
            'tol' : [1e-2, 5e-3, 1e-3],
            'ts' : ['midpoint', 'bdf2'],
            }

    elif parameter_set == 'cubeoid-timestep-newton-convergence':
        args_dict = {
            'driver' : ["./llg_driver/llg_driver"],
            'dt' : [1e-4, 1e-5, 1e-6, 1e-7, 1e-8],
            'tmax' : [1e-8],
            'tol' : [0.0],
            'ref' : [3],
            'ts' : ['bdf2'],
            'initm' : ['z'],
            'happ' : ['minus_z'],
            'mesh' : ['ut_cubeoid', 'sq_cubeoid', 'st_cubeoid'],
            'solver' : ['superlu'],
            }

    elif parameter_set == 'oscillating_fields':
        args_dict = {
            'driver' : ["./llg_driver/llg_driver"],
            'dt' : [1e-4],
            'tmax' : [200],
            'tol' : [1e-3, 5e-4, 1e-4, 1e-5],
            'ref' : [2,3,4],
            'ts' : ['bdf2', 'midpoint'],
            'initm' : ['z'],
            'happ' : ['z_oscillating_p20'],
            'mesh' : ['sq_square'],
            'mag-params' :["simple-llg-max-damped"]
            }


    elif parameter_set == "prec_fast":
        args_dict = {
            'driver' : ["./llg_driver/llg_driver"],
            'dt' : [0.01],
            'ref' : [2],
            'tmax' : [1.0],
            'ts' : ["bdf2"],
            'initm' : ['smoothly_varying_50'],
            'happ' : ['x', 'y', 'z'],
            'mesh' : ['sq_square', 'ut_square'],
            'solver' : ['gmres'],
            'preconditioner': ['blockllg-uppertriangular-blockantidiagonal-xy-exact',
                               'blockllg-uppertriangular-blockantidiagonal-xz-exact',
                               'blockllg-uppertriangular-blockantidiagonal-yz-exact',
                               ]

            }

    elif parameter_set == "prec":
        args_dict = {
            'driver' : ["./llg_driver/llg_driver"],
            'dt' : [0.1, 0.05, 0.01, 0.001],
            'ref' : [1, 2, 3, 4],
            'tmax' : [1.0],
            'ts' : ["bdf2"],
            'initm' : ['smoothly_varying_50'],
            'happ' : ['x', 'y', 'z', 'all_directions'],
            'mesh' : ['sq_square', 'ut_square', 'ut_sphere'],
            'solver' : ['gmres'],
            'preconditioner': ['blockllg-uppertriangular-blockexact-xy-exact',
                               'blockllg-uppertriangular-blockexact-xz-exact',
                               'blockllg-uppertriangular-blockexact-yz-exact',
                               ]
            }

    elif parameter_set == "adaptive-midpoint-conservation":
        args_dict = {
            'driver' : ["./semi_implicit_mm_driver/semi_implicit_mm_driver"],
            'dt' : [1e-5],
            'tol' : [1e-2, 1e-3, 1e-4],
            'ref' : [3, 4, 5],
            'tmax' : [10.0],
            'ts' : ["bdf2", "midpoint"],
            'initm' : ['smoothly_varying_50'],
            'happ' : ['minus_z'],
            'mesh' : ['ut_square'],
            'renorm_m' : [0]
            }

    elif parameter_set == "fixed-step-midpoint-conservation":
        args_dict = {
            'driver' : ["./semi_implicit_mm_driver/semi_implicit_mm_driver"],
            'dt' : [1e-1, 1e-2, 1e-3, 5e-3],
            'ref' : [3],
            'tmax' : [5.0],
            'ts' : ["bdf2", "midpoint"],
            'initm' : ['z'],
            'happ' : ['minus_z'],
            'mesh' : ['ut_sphere'],
            'renorm_m' : [1]
            }

    elif parameter_set == "zero-damping":
        args_dict = {
            'driver' : ["./semi_implicit_mm_driver/semi_implicit_mm_driver"],
            'dt' : [1e-1, 1e-2, 1e-3, 5e-3],
            'ref' : [3],
            'tmax' : [5.0],
            'ts' : ["bdf2", "midpoint"],
            'mesh' : ['sq_square'],
            'renorm_m' : [1],
            'dampc' : [0]
            }

    elif parameter_set == "adaptive-midpoint":
        args_dict = {
            'driver' : ["./llg_driver/llg_driver"],
            'tol' : [1e-2, 1e-3, 1e-4],
            'ref' : [1],
            'tmax' : [4.0],
            'ts' : ["bdf2", "midpoint"],
            'mesh' : ['sq_square'],
            'renorm_m' : [1],
            'dampc' : [0.5],
            'resi' : ["ll"]
            }

    elif parameter_set == "ode-test":
         args_dict = {
            'driver' : ["./ode_problem/ode_problem"],
            'tol' : [1e-2, 1e-3],
            'tmax' : [10.0],
            'ts' : ["bdf2", "midpoint"],
            'mp-pred' : ["edbdf3", "rk4"]
            }

    elif parameter_set == "multi-domain-failures":
         args_dict = {
            'driver' : ["./semi_implicit_mm_driver/semi_implicit_mm_driver"],
            'dt' : [1e-4, 1e-5],
            'tmax' : [1.0],
            'ts' : ["bdf2", "midpoint"],
            'mesh' : ['multi_ut_square', 'multi_sq_square'],
            'renorm_m' : [1],
            'ref' : [3],
            'doc-interval' : ["all"]
            }

    elif parameter_set == "multi-domain-squares":
         args_dict = {
            'driver' : ["./semi_implicit_mm_driver/semi_implicit_mm_driver"],
            'tmax' : [1e-10],
            'ts' : ["midpoint"],
            'mesh' : ['multi_sq_square', 'sq_square'],
            'ref' : [1],
            'happ' : ['zero'],
            }

    elif parameter_set == "multi-domain-spheres-energy-test":
         args_dict = {
            'driver' : ["./semi_implicit_mm_driver/semi_implicit_mm_driver"],
            'tmax' : [1e-10],
            'ts' : ["midpoint"],
            'mesh' : ['multi_ut_sphere', 'ut_sphere'],
            'ref' : [3],
            'happ' : ['zero'],
            }

    elif parameter_set == "compare-implicitness-semi":
         args_dict = {
            'driver' : ["./semi_implicit_mm_driver/semi_implicit_mm_driver"],
            'tmax' : [20],
            'ts' : ["bdf2"],
            'mesh' : ['many_ut_square'],
            'ref' : [3, 4],
            'tol' : [1e-1, 1e-3, 1e-5],
            'dt' : [1e-6],
            }

    elif parameter_set == "compare-implicitness-implicit":
         args_dict = {
            'driver' : ["./llg_driver/llg_driver"],
            'tmax' : [20],
            'ts' : ["bdf2"],
            'mesh' : ['many_ut_square'],
            'ref' : [3, 4],
            'tol' : [1e-1, 1e-3, 1e-5],
            'dt' : [1e-6],
            'implicit-ms' : [True],
            'solver' : ['som-gmres'],
            'preconditioner' : ['som-main-exact']
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

    # Make sure micromag library is up to date
    library_folder = os.path.abspath("../")
    print("Building and installing libraries from", library_folder)
    subp.check_call(['make', '--silent', '--keep-going',
                     'LIBTOOLFLAGS=--silent'], cwd=library_folder)
    subp.check_call(['make', 'install', '--silent', '--keep-going',
                     'LIBTOOLFLAGS=--silent'], cwd=library_folder)

    # Make sure the driver binaries are up to date
    driver_folders = [os.path.dirname(d) for d in args_dict['driver']]
    map(build, driver_folders)

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

    parser.add_argument('--clean', action='store_true',
                        help='clean up old results from the target folder')

    # parser.add_argument('-ncores', '-j', dest='ncores',
    #                     help='Set number of cores to use.')

    args = parser.parse_args()

    # Main function
    # ============================================================


    # Do parameter sweep
    if args.parameters is not None:
        standard_sweep(args.parameters, args.clean, args.debug_mode)

    # Or just dump some Jacobians
    elif args.j_parameter_set is not None:
        print("Running Jacobian generation parameter sweep",
              "with parameter set", args.j_parameter_set)
        milan_jacobians(args.j_parameter_set,args.debug_mode)

    return 0


# If this script is run from a shell then run main() and return the result.
if __name__ == "__main__":
    sys.exit(main())
