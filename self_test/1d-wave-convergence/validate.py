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
        # Problem specification
        '-driver' : 'llg',
        '-ms-method' : 'disabled',
        '-mesh' : 'sq_line_periodic',
        '-initial-m' : 'periodic_exact',
        '-initial-is-exact' : 1,
        '-h-app' : 'zero',
        '-damping' : [0.9, 0.1, 0.01, 0.001, 0],
        '-tmax' : 0.1,
        '-wave-solution-c' : 1/12, # as used by Jeong et. al.

        # convergence test: one step and link dt to spatial refinement
        '-max-steps' : 1,
        '-convergence-test' : 1,
        '-doc-interval' : 0,
        '-doc-exact' : 1,


        # Integration/calculation details
        '-ts' : ["imr", "tr", "bdf2"],
        '-ref' : [2, 3, 4, 5, 6, 7, 8],
        '-newton-tol' : 1e-12,
        '-renormalise' : [0],
        '-quadrature' : ['lnodal', 'gauss'],
        }

    # Where it's going to end up
    base_outdir = os.path.abspath(pjoin(os.path.dirname(__file__), "Validation"))

    # Run
    err_codes, outdirs = mm.run_sweep(argdicts, base_outdir,
                                      parallel_sweep=args.parallel)

    # Get data
    datasets = list(map(mm.parse_run, outdirs))

    # Check things all ran
    ran = all((e == 0 for e in err_codes))

    convergence_test_datasets = mm.split_to_comparable_groups(datasets, '-ref')

    def rate(datasets):
        """Expected convergence rate for a given timestepper.
        """
        # Check that all have the same ts (they should since this is a
        # convergence test!)
        assert all((d['-ts'] == datasets[0]['-ts'] for d in datasets))

        # bdf2 convergence is pretty bad, choose a lower convergence rate for it
        if datasets[0]['-ts'] == 'bdf2':
            return 1.75
        else:
            return 2

    # Check the convergence rates, seem to be just over 2, not sure why
    # since it should only be 2. Something to do with using an exact
    # solution?
    t1 = all([tests.check_convergence(d, 2.2, tol=0.2)
              for d in convergence_test_datasets])

    if ran and t1:
        return 0
    else:
        return 1


if __name__ == "__main__":
    sys.exit(main())
