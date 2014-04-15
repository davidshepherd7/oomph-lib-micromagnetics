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
    """Check that all the preconditioners run without crashing (doesn't check that they work robustly because they don't yet!)
    """

    # Look for parallel in args
    parser = argparse.ArgumentParser()
    parser.add_argument('--parallel', action = "store_true")
    args = parser.parse_args()

    llg_precs = ["exact", "block"]
    llg_sub_precs = ["exact", "block"]

    # With magnetostatics
    argdicts_ms = {
        "-driver" : 'llg',
        "-dt" : 0.1,
        "-solver" : "som-gmres",
        "-prec" : "som-main-block",
        "-llg-prec" : llg_precs,
        "-llg-sub-prec" : llg_sub_precs,
        }

    # Without magnetostatics
    argdicts = {
        "-driver" : 'llg',
        "-dt" : 0.1,
        "-ms-method" : "disabled",
        "-solver" : "gmres",
        "-prec" : "dummy",
        "-llg-prec" : llg_precs,
        "-llg-sub-prec" : llg_sub_precs,
        }

    # Slightly redundant to run with multiple llg-sub-prec args even when
    # using exact llg-prec (llg-sub-prec is ignored in this case). But it's
    # so fast that it really doesn't matter.

    # Where it's going to end up
    base_outdir = os.path.abspath(pjoin(os.path.dirname(__file__), "Validation"))

    # Run
    err_codes_ms, _ = mm.run_sweep(argdicts_ms, base_outdir,
                                   parallel_sweep=args.parallel)
    err_codes, _ = mm.run_sweep(argdicts, base_outdir,
                                parallel_sweep=args.parallel)

    # Check things ran ok
    ran = all((e == 0 for e in it.chain(err_codes_ms, err_codes)))

    # No analysis of output: not checking that here

    if ran:
        return 0
    else:
        return 1


if __name__ == "__main__":
    sys.exit(main())
