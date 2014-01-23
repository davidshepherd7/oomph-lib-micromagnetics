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


def build_driver(folder):
    print("Building in", folder)
    subp.check_call(['make', '--silent', '--keep-going',
                    'LIBTOOLFLAGS=--silent'], cwd=folder)


def milan_jacobians(parameter_set, serial_mode=False):

    # Presumably we will always be using the implicit driver...


    # Construct lists of args
    if parameter_set == "initial":
        args_dict = {
            '-binary' : "./llg_driver/llg_driver",
            '-dt' : [0.1, 0.05, 0.01, 0.001],
            '-tmax' : [0.001],
            '-tol' : [0.0],
            '-ref' : [1, 2, 3, 4, 5],
            '-ts' : ["bdf2"],
            '-initm' : ['smoothly_varying'],
            '-happ' : ['x', 'y', 'z'],
            '-mesh' : ['sq_square', 'ut_square'],
            '-output-jac' : ['at_end']
            }
        # more parameter sets go here
    else:
        raise NotImplementedError

    # Fill in some function args that should always be the same.
    fun = par(execute_oomph_driver,
              output_root='../experiments/jacobian_sweeps')

    # Run the parameter sweep!
    mm.parallel_parameter_sweep(fun, args_dict, serial_mode)

    return 0

def standard_sweep(parameter_set, cleanup, serial_mode=False):

    if parameter_set == 'nmag_cubeoid':
        args_dict = {
            '-binary' : ["./semi_implicit_mm_driver/semi_implicit_mm_driver"],
            '-dt' : [1e-4],
            '-tmax' : [60.0],
            '-tol' : [1e-3, 1e-4, 1e-5],
            '-ref' : [4],
            '-ts' : ['midpoint', 'bdf2'],
            '-initm' : ['xz'],
            '-happ' : ['zero'],
            '-mesh' : ['sq_cubeoid'],
            '-solver' : ['gmres'],
            '-prec' : ['amg']
            }

    elif parameter_set == 'const_dt_nmag_cubeoid':
        args_dict = {
            '-binary' : ["./semi_implicit_mm_driver/semi_implicit_mm_driver"],
            '-dt' : [1e-2, 5e-3],
            '-tmax' : [20.0],
            '-tol' : [0.0],
            '-ref' : [2, 3, 4],
            '-ts' : ['midpoint', 'bdf2'],
            '-initm' : ['xz'],
            '-happ' : ['zero'],
            '-mesh' : ['sq_cubeoid'],
            '-solver' : ['gmres'],
            '-prec' : ['amg']
            }

    elif parameter_set == 'nmag_cubeoid_llg_prec':
        args_dict = {
            '-binary' : ["./semi_implicit_mm_driver/semi_implicit_mm_driver"],
            '-dt' : [5e-4],
            '-tmax' : [20.0],
            '-tol' : [1e-4],
            '-ref' : [4],
            '-ts' : ['bdf2'],
            '-initm' : ['xz'],
            '-happ' : ['zero'],
            '-mag-params' :["simple-llg"],
            '-mesh' : ['sq_cubeoid'],
            '-solver' : ['gmres'],
            '-prec' : ['blockllg-uppertriangular-blockexact-xy-exact',
                                'blockllg-uppertriangular-blockexact-xz-exact',
                                'blockllg-uppertriangular-blockexact-yz-exact',
                                ]
            }

    elif parameter_set == 'script_test':
        args_dict = {
            '-driver' : ["ode"],
            '-dt' : [1e-4],
            '-tmax' : [1.0],
            '-tol' : [1e-3],
            '-ref' : [2],
            }

    elif parameter_set == 'unsteady_heat_midpoint_vs_bdf2':
        args_dict = {
            '-binary' : ["./unsteady_heat_driver/unsteady_heat_driver"],
            '-dt' : [1e-4],
            '-tmax' : [10.0],
            '-tol' : [1e-2, 5e-3, 1e-3],
            '-ts' : ['midpoint', 'bdf2'],
            }

    elif parameter_set == 'cubeoid-timestep-newton-convergence':
        args_dict = {
            '-binary' : ["./llg_driver/llg_driver"],
            '-dt' : [1e-4, 1e-5, 1e-6, 1e-7, 1e-8],
            '-tmax' : [1e-8],
            '-tol' : [0.0],
            '-ref' : [3],
            '-ts' : ['bdf2'],
            '-initm' : ['z'],
            '-happ' : ['minus_z'],
            '-mesh' : ['ut_cubeoid', 'sq_cubeoid', 'st_cubeoid'],
            '-solver' : ['superlu'],
            }

    elif parameter_set == 'oscillating_fields':
        args_dict = {
            '-binary' : ["./llg_driver/llg_driver"],
            '-dt' : [1e-4],
            '-tmax' : [200],
            '-tol' : [1e-3, 5e-4, 1e-4, 1e-5],
            '-ref' : [2,3,4],
            '-ts' : ['bdf2', 'midpoint'],
            '-initm' : ['z'],
            '-happ' : ['z_oscillating_p20'],
            '-mesh' : ['sq_square'],
            '-mag-params' :["simple-llg-max-damped"]
            }


    elif parameter_set == "prec_fast":
        args_dict = {
            '-binary' : ["./llg_driver/llg_driver"],
            '-dt' : [0.01],
            '-ref' : [2],
            '-tmax' : [1.0],
            '-ts' : ["bdf2"],
            '-initm' : ['smoothly_varying_50'],
            '-happ' : ['x', 'y', 'z'],
            '-mesh' : ['sq_square', 'ut_square'],
            '-solver' : ['gmres'],
            '-prec': ['blockllg-uppertriangular-blockantidiagonal-xy-exact',
                     'blockllg-uppertriangular-blockantidiagonal-xz-exact',
                     'blockllg-uppertriangular-blockantidiagonal-yz-exact',
                     ]

            }

    elif parameter_set == "prec":
        args_dict = {
            '-binary' : ["./llg_driver/llg_driver"],
            '-dt' : [0.1, 0.05, 0.01, 0.001],
            '-ref' : [1, 2, 3, 4],
            '-tmax' : [1.0],
            '-ts' : ["bdf2"],
            '-initm' : ['smoothly_varying_50'],
            '-happ' : ['x', 'y', 'z', 'all_directions'],
            '-mesh' : ['sq_square', 'ut_square', 'ut_sphere'],
            '-solver' : ['gmres'],
            '-prec': ['blockllg-uppertriangular-blockexact-xy-exact',
                     'blockllg-uppertriangular-blockexact-xz-exact',
                     'blockllg-uppertriangular-blockexact-yz-exact',
                     ]
            }

    elif parameter_set == "adaptive-midpoint-conservation":
        args_dict = {
            '-binary' : ["./semi_implicit_mm_driver/semi_implicit_mm_driver"],
            '-dt' : [1e-5],
            '-tol' : [1e-2, 1e-3, 1e-4],
            '-ref' : [3, 4, 5],
            '-tmax' : [10.0],
            '-ts' : ["bdf2", "midpoint"],
            '-initm' : ['smoothly_varying_50'],
            '-happ' : ['minus_z'],
            '-mesh' : ['ut_square'],
            '-renorm_m' : [0]
            }

    elif parameter_set == "fixed-step-midpoint-conservation":
        args_dict = {
            '-binary' : ["./semi_implicit_mm_driver/semi_implicit_mm_driver"],
            '-dt' : [1e-1, 1e-2, 1e-3, 5e-3],
            '-ref' : [3],
            '-tmax' : [5.0],
            '-ts' : ["bdf2", "midpoint"],
            '-initm' : ['z'],
            '-happ' : ['minus_z'],
            '-mesh' : ['ut_sphere'],
            '-renorm_m' : [1]
            }

    elif parameter_set == "zero-damping":
        args_dict = {
            '-binary' : ["./semi_implicit_mm_driver/semi_implicit_mm_driver"],
            '-dt' : [1e-1, 1e-2, 1e-3, 5e-3],
            '-ref' : [3],
            '-tmax' : [5.0],
            '-ts' : ["bdf2", "midpoint"],
            '-mesh' : ['sq_square'],
            '-renorm_m' : [1],
            '-dampc' : [0]
            }

    elif parameter_set == "adaptive-midpoint":
        args_dict = {
            '-binary' : ["./llg_driver/llg_driver"],
            '-tol' : [1e-2, 1e-3, 1e-4],
            '-ref' : [1],
            '-tmax' : [4.0],
            '-ts' : ["bdf2", "midpoint"],
            '-mesh' : ['sq_square'],
            '-renorm_m' : [1],
            '-dampc' : [0.5],
            '-resi' : ["ll"]
            }

    elif parameter_set == "ode-test":
         args_dict = {
            '-binary' : ["./ode_problem/ode_problem"],
            '-tol' : [1e-2, 1e-3],
            '-tmax' : [10.0],
            '-ts' : ["bdf2", "midpoint"],
            '-mp-pred' : ["edbdf3", "rk4"]
            }

    elif parameter_set == "multi-domain-failures":
         args_dict = {
            '-binary' : ["./semi_implicit_mm_driver/semi_implicit_mm_driver"],
            '-dt' : [1e-4, 1e-5],
            '-tmax' : [1.0],
            '-ts' : ["bdf2", "midpoint"],
            '-mesh' : ['multi_ut_square', 'multi_sq_square'],
            '-renorm_m' : [1],
            '-ref' : [3],
            '-doc-interval' : ["all"]
            }

    elif parameter_set == "multi-domain-squares":
         args_dict = {
            '-binary' : ["./semi_implicit_mm_driver/semi_implicit_mm_driver"],
            '-tmax' : [1e-10],
            '-ts' : ["midpoint"],
            '-mesh' : ['multi_sq_square', 'sq_square'],
            '-ref' : [1],
            '-happ' : ['zero'],
            }

    elif parameter_set == "multi-domain-spheres-energy-test":
         args_dict = {
            '-binary' : ["./semi_implicit_mm_driver/semi_implicit_mm_driver"],
            '-tmax' : [1e-10],
            '-ts' : ["midpoint"],
            '-mesh' : ['multi_ut_sphere', 'ut_sphere'],
            '-ref' : [3],
            '-happ' : ['zero'],
            }

    elif parameter_set == "compare-implicitness-semi":
         args_dict = {
            '-binary' : ["./semi_implicit_mm_driver/semi_implicit_mm_driver"],
            '-tmax' : [20],
            '-ts' : ["bdf2"],
            '-mesh' : ['many_ut_square'],
            '-ref' : [3, 4],
            '-tol' : [1e-1, 1e-3, 1e-5],
            '-dt' : [1e-6],
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

    # Make sure micromag library is up to date
    library_folder = os.path.abspath("../")
    print("Building and installing libraries from", library_folder)
    subp.check_call(['make', '--silent', '--keep-going',
                     'LIBTOOLFLAGS=--silent'], cwd=library_folder)
    subp.check_call(['make', 'install', '--silent', '--keep-going',
                     'LIBTOOLFLAGS=--silent'], cwd=library_folder)

    # Make sure the driver binaries are up to date (if they aren't just the
    # default)
    binaries = args_dict.get('-binary')
    if binaries is not None:
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
