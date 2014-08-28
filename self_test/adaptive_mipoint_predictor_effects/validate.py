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


def main():
    """Check that using adaptivity doesn't modify results if timestep changes
    are ignored.
    """

    base_outdir = os.path.abspath(pjoin(os.path.dirname(__file__), "Validation"))

    # What to run
    base_argdict = {
        "-driver" : ["llg"],
        "-ts" : "midpoint-bdf",
        "-ms-method" : ["implicit", "decoupled", "disabled"],
        "-solver" : "gmres",
        "-matrix-type" : "som",
        "-prec" : "som-main-exact",
        "-tmax" : 15,
        '-dt-initial' : 0.1,
        '-newton-tol' : 1e-12,
        }

    fa_argdict = copy.deepcopy(base_argdict)
    fa_argdict.update({"-tol" : 1e-3, '-dummy-adaptivity' : 1})
    fa_err_codes, fa_outdirs = mm.run_sweep(fa_argdict, pjoin(base_outdir, "fa"))


    cs_argdict = copy.deepcopy(base_argdict)
    cs_argdict.update({"-tol" : 0})
    cs_err_codes, cs_outdirs = mm.run_sweep(cs_argdict, pjoin(base_outdir, "cs"))

    ran = all([e == 0 for e in fa_err_codes + cs_err_codes])

    with open(os.devnull, 'w') as null:
        ok = all([tests.check_solns_match(a, b, outstream=null,
                                          details_stream=null)
                  for a, b in zip(fa_outdirs, cs_outdirs)])

    if ran and ok:
        return 0
    else:
        return 1


if __name__ == "__main__":
    sys.exit(main())
