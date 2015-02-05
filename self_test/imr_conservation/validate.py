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
    parser.add_argument('--short', action = "store_true")
    args = parser.parse_args()


    # What to run
    argdicts_1d = {
         # Problem specification
        '-driver' : 'llg',
        '-ms-method' : 'disabled',
        '-mesh' : 'sq_line_periodic',
        '-initial-m' : 'periodic_exact',
        '-initial-is-exact' : 1,
        '-h-app' : 'zero',
        '-damping' : [0.9, 0.01, 0],
        '-tmax' : 1.0,
        '-wave-solution-c' : 1/12, # Jeong et. al's value


        # Integration/calculation details
        '-ts' : ["imr"],
        '-ref' : [1, 4],
        '-dt' : [0.1, 0.01, 0.005],
        # ref 1 is way to little and dt 0.1 is way to big but test
        # conservation properties anyway (should conserve).

        '-newton-tol' : 1e-12,
        '-renormalise' : "never",
        '-quadrature' : ['lnodal'],
        }

    argdicts_2d = {
         # Problem specification
        '-driver' : 'llg',
        '-ms-method' : 'disabled',
        '-mesh' : 'sq_square_periodic',
        '-initial-m' : 'periodic_exact',
        '-initial-is-exact' : 1,
        '-h-app' : 'zero',
        '-damping' : [0.5], #only one case because it's slower
        '-tmax' : 1.0,
        '-wave-solution-c' : 1/12, # Jeong et. al's value

        # Integration/calculation details
        '-ts' : ["imr"],
        '-ref' : [2, 5], # need to test high refinement to check that
                         # rescaling is working
        '-dt' : [0.1],
        '-newton-tol' : 1e-12,
        '-renormalise' : "never",
        '-quadrature' : ['lnodal'],
        }


    # Seems like exact solution is weird, so check another case as well
    argdicts_non_exact = {
         # Problem specification
        '-driver' : 'llg',
        '-ms-method' : 'disabled',
        '-mesh' : ['sq_square'], # 'st_square'],
        '-initial-m' : 'smoothly_varying_5',
        '-h-app' : 'zero',
        '-damping' : [0.5],
        '-tmax' : 1.0,

        # Integration/calculation details
        '-ts' : ["imr"],
        '-ref' : [2],
        '-dt' : [0.1],
        '-newton-tol' : 1e-12,
        '-renormalise' : "never",
        '-quadrature' : ['lnodal'],
        }


    # Where it's going to end up
    base_outdir = os.path.abspath(pjoin(os.path.dirname(__file__), "Validation"))

    # Run
    if not args.short:
        err_codes_1d, outdirs_1d = mm.run_sweep(argdicts_1d, base_outdir,
                                                serial_mode=not args.parallel)
        err_codes_2d, outdirs_2d = mm.run_sweep(argdicts_2d, base_outdir + "_2d",
                                                serial_mode=not args.parallel)
    else:
        err_codes_1d = []
        outdirs_1d = []
        err_codes_2d = []
        outdirs_2d = []

    err_codes_non_exact, outdirs_non_exact = mm.run_sweep(argdicts_non_exact,
                                                          base_outdir+"_non_exact",
                                                          serial_mode=not args.parallel)

    err_codes = err_codes_1d + err_codes_2d + err_codes_non_exact
    outdirs = outdirs_1d + outdirs_2d + outdirs_non_exact


    # Get data
    datasets = list(map(mm.parse_run, outdirs))

    # Check things
    ran = all((e == 0 for e in err_codes))

    t1 = all([tests.check_m_length(d, tol=1e-10) for d in datasets])

    if ran and t1:
        return 0
    else:
        return 1


if __name__ == "__main__":
    sys.exit(main())
