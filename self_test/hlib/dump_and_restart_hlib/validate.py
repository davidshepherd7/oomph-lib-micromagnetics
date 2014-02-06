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
import oomphpy.micromagnetics as mm
import oomphpy.tests as tests


def main():

    argdicts = {
        "-driver" : ["llg"],
        '-dump' : [1],
        '-dt' : [0.1],
        '-max-steps' : [30],
        '-tmax' : [999],
        '-mesh' : ['ut_cubeoid', 'st_cubeoid'],
        '-decoupled-ms' : [True, False],
        '-doc-interval' : [0],
        '-hlib-bem' : [1],
        '-solver' : ['som-gmres'],
        '-prec' : ['som-main-exact'],
        }

    # Where it's going to end up
    base_outdir = os.path.abspath(pjoin(os.path.dirname(__file__), "Validation",
                                        "initial"))
    base_restart_outdir = os.path.abspath(pjoin(os.path.dirname(__file__), "Validation",
                                        "restart"))

    # Run
    err_codes, outdirs = mm.run_sweep(argdicts, base_outdir)


    varying_args = []
    for k, v in argdicts.items():
        if len(v) > 1:
            varying_args.append(k)


    # Run again from restart, write loop rather than using run_sweep so
    # that we can handle restarts
    restart_err_codes = []
    restart_outdirs = []
    # ??ds massive hack here: I'm hoping that the order of outdirs is the
    # same as the order of mm.product_of_argdict(argdicts)...
    for restart_argdict, old_outdir in zip(mm.product_of_argdict(argdicts),
                                            outdirs):
        restart_argdict['-restart'] = os.path.abspath(pjoin(old_outdir, "dump20.dat"))

        restart_argdict['-dump'] = 0

        err, outdir = mm._run(restart_argdict, base_restart_outdir,
                              varying_args)

        restart_err_codes.append(err)
        restart_outdirs.append(outdir)


    # Check things
    ran = all((e == 0 for e in it.chain(err_codes, restart_err_codes)))

    t1 = all([tests.check_restarted_in_middle(rdir, 20) for rdir in restart_outdirs])

    with open(os.devnull, 'w') as null:
        fpdiff_args = {'details_stream' : null,
                       # outstream=null,
                       'relative_error' : 0.1,
                       'small' : 1e-14,
                       }

        t2 = all([tests.check_solns_match(rdir, odir, **fpdiff_args)
                  for rdir, odir in zip(restart_outdirs, outdirs)])

    if ran and t1 and t2:
        return 0
    else:
        return 1


if __name__ == "__main__":
    sys.exit(main())
