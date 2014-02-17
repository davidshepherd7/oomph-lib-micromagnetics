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

    # Test without magnetostatics by comparison with Mallinson solution
    noms_argdicts = {
        "-driver" : ["ll"],
        "-dt": [0.01],
        "-scale": [10],
        "-ms-method" : "disabled",
        "-ts" : ["cay-euler", "cay-rk2"],
        "-tmax" : [3],
        }

    base_outdir = os.path.abspath(pjoin(os.path.dirname(__file__), "Validation"))

    # Run
    err_codes, outdirs, _ = mm.run_sweep(noms_argdicts, base_outdir)

    # Get data
    datasets = list(map(mm.parse_run, outdirs))

    # Check things
    ran = all((e == 0 for e in err_codes))
    t1 = all([tests.check_error_norm(d, 1e-5) for d in datasets])
    t2 = all([tests.check_m_length(d, 1e-14) for d in datasets])

    if ran and t1 and t2:
        return 0
    else:
        return 1



if __name__ == "__main__":
    sys.exit(main())
