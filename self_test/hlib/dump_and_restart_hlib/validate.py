#!/usr/bin/env python3

# Python 2/3 compatability
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import sys
import argparse
import os
import os.path

import itertools as it

from os.path import join as pjoin

# Make sure *this* versions oomphpy is in the path (before any other
# versions in other places)
sys.path.insert(1, pjoin(os.path.dirname(__file__), "../../../etc"))
import oomphpy
import oomphpy.utils as utils
import oomphpy.micromagnetics as mm
import oomphpy.tests as tests


def restarted_run(restart_argdict, varying_args, base_restart_outdir, base_outdir):

    # Construct original outdir name
    varying_arg_names = [str(restart_argdict[k]) for k in sorted(varying_args)]
    old_outdir = pjoin(base_outdir, "results_" + "_".join(varying_arg_names))

    restart_argdict['-restart'] = os.path.abspath(pjoin(old_outdir, "dump2.dat"))
    restart_argdict['-dump'] = 0
    return mm._run(restart_argdict, base_restart_outdir, varying_args)


def main():

    # Parse arguments
    # ============================================================

    parser = argparse.ArgumentParser(description=main.__doc__,
        # Don't mess up my formating in the help message
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--parallel', '-p', action = "store_true")
    parser.add_argument('--slow', '-s', action = "store_true")
    args = parser.parse_args()

    if args.slow:
        max_step = 10
        meshes = ['ut_cubeoid', 'st_cubeoid', 'ut_sphere']
    else:
        max_step = 4
        meshes = ['ut_cubeoid']

    argdicts = {
        "-driver" : ["llg"],
        '-dump' : [1],
        '-dt' : [0.1],
        '-max-steps' : [max_step],
        '-tmax' : [999],
        '-mesh' : meshes,
        '-ms-method' : ["decoupled", "disabled"],
        '-doc-interval' : [0],
        '-hlib-bem' : [1],
        '-solver' : ['som'],
        "-matrix-type" : "som",
        '-prec' : ['som-main-ilu-0'],
        '-phi1-singularity-method' : 'pin_bulk',
        }

    # Where it's going to end up
    base_outdir = os.path.abspath(pjoin(os.path.dirname(__file__), "Validation",
                                        "initial"))
    base_restart_outdir = os.path.abspath(pjoin(os.path.dirname(__file__), "Validation",
                                        "restart"))

    # Run
    err_codes, outdirs = mm.run_sweep(argdicts, base_outdir, args.parallel)


    varying_args = []
    for k, v in argdicts.items():
        if len(v) > 1:
            varying_args.append(k)


    # Run again from restart
    restart_err_codes, restart_outdirs = \
        utils.unzip( utils.parallel_map(restarted_run, mm.product_of_argdict(argdicts),
                                it.repeat(varying_args),
                                it.repeat(base_restart_outdir),
                                it.repeat(base_outdir),
                                serial_mode=not args.parallel) )

    # Check things
    ran = all((e == 0 for e in it.chain(err_codes, restart_err_codes)))

    t1 = all([tests.check_restarted_in_middle(rdir, 2) for rdir in restart_outdirs])

    with open(os.devnull, 'w') as null:
        fpdiff_args = {'details_stream' : null,
                       'outstream' : null,
                       'relative_error' : 0.1,
                       'small' : 1e-8,
                       }

        t2 = all([tests.check_solns_match(rdir, odir, **fpdiff_args)
                  for rdir, odir in zip(restart_outdirs, outdirs)])

    if ran and t1 and t2:
        return 0
    else:
        return 1


if __name__ == "__main__":
    sys.exit(main())
