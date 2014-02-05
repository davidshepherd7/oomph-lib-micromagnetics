#!/usr/bin/env python3

# Python 2/3 compatability
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import sys
import argparse
import os
import os.path
import copy

import itertools as it

from os.path import join as pjoin

# Make sure *this* versions oomphpy is in the path (before any other
# versions in other places)
sys.path.insert(1, pjoin(os.path.dirname(__file__), "../../etc"))
import oomphpy
import oomphpy.micromagnetics as mm
import oomphpy.tests as tests

from fpdiff import fpdiff


def check_restarted_in_middle(restart_outdir):
    ok = (not os.path.isfile(pjoin(restart_outdir, "soln19.dat"))) \
      and (os.path.isfile(pjoin(restart_outdir, "soln20.dat")))

    if ok:
        mm.okprint("Restarted at correct point in", restart_outdir)
    else:
        mm.badprint("Failed or restarted at wrong point in", restart_outdir)

    return ok

def check_restart_solns_correct(restart_outdir, outdir):

    matches = []
    for fname in os.listdir(restart_outdir):
        if ("soln" in fname) and (".dat" in fname):
            matches.append(fpdiff(pjoin(restart_outdir, fname),
                                  pjoin(outdir, fname),
                                  0.1,
                                  1e-14))

    if len(matches) == 0:
        mm.badprint("No files in", restart_outdir)
        return False

    ok = all(matches)
    if ok:
        mm.okprint("Files match in", restart_outdir, outdir)
    else:
        mm.badprint("Files don't match in", restart_outdir, outdir)

    return ok


def main():

    # What to run
    argdicts = {
        "-driver" : ["llg"],
        '-dump' : [1],
        '-dt' : [0.1],
        '-max-steps' : [30],
        '-tmax' : [999],
        '-mesh' : ['sq_square'],# , 'ut_square'],
        '-decoupled-ms' : [False],
        '-disable-ms' : [True],
        '-doc-interval' : [0],
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
    t2 = all([check_restart_solns_correct(rdir, odir)
              for rdir, odir in zip(restart_outdirs, outdirs)])
    t1 = all([check_restarted_in_middle(rdir) for rdir in restart_outdirs])

    if ran and t1 and t2:
        return 0
    else:
        return 1


if __name__ == "__main__":
    sys.exit(main())
