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
sys.path.insert(1, pjoin(os.path.dirname(__file__), "../../etc"))
import oomphpy
import oomphpy.micromagnetics as mm
import oomphpy.tests as tests


def main():

    # Look for parallel in args
    parser = argparse.ArgumentParser()
    parser.add_argument('--parallel', action = "store_true")
    args = parser.parse_args()


    # What to run
    argdicts = {
        "-driver" : 'llg',

        '-mesh' : 'ut_sphere',
        '-tmax' : 10,
        '-damping' : 1,

        '-ts' : 'imr',
        '-dt' : 0.1,

        '-ms-method' : 'implicit',
        '-quadrature' : 'lnodal',
        '-newton-tol' : 1e-12,
        "-solver" : "gmres",
        "-matrix-type" : "som",
        "-prec" : "som-main-exact",
        }

    # Where it's going to end up
    base_outdir = os.path.abspath(pjoin(os.path.dirname(__file__), "Validation"))

    # Run
    err_codes, outdirs = mm.run_sweep(argdicts, base_outdir,
                                      serial_mode=not args.parallel)

    # Get data
    datasets = list(map(mm.parse_run, outdirs))

    # Check things
    ran = all((e == 0 for e in err_codes))
    t1 = all([tests.check_m_length(d, 1e-7) for d in datasets])
    t2 = all([tests.check_error_norm(d, 1e-2) for d in datasets])

    # ??ds not convinced I should let the tols be this loose

    if ran and t1:
        return 0
    else:
        return 1


if __name__ == "__main__":
    sys.exit(main())
