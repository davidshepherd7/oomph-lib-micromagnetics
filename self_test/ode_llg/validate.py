#!/usr/bin/env python3

# Python 2/3 compatability
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import sys
import argparse
import os
import os.path

from os.path import join as pjoin

# Make sure *this* versions oomphpy is in the path (before any other
# versions in other places)
sys.path.insert(1, pjoin(os.path.dirname(__file__), "../../etc"))
import oomphpy
import oomphpy.micromagnetics as mm


def small(f):
    return abs(f) < 1e-4

def main():

    # ??ds ad "llg"?
    for exact in ["ll"]:
        for ts in ["rk2", "midpoint-bdf"]:

            # parameters
            valdir = pjoin(".", "Validation", exact+"_"+ts)
            args = ['llgode', '-exact', exact,
                    '-ts', ts,
                    '-dt', 0.01,
                    '-outdir', valdir]
            # Run
            mm.cleandir(valdir)
            err_code = mm.run_driver(args, valdir)

            assert(err_code == 0)

            # Check all errors are small
            data = mm.parse_trace_file(pjoin(valdir, "trace"))


            max_err = max(map(abs, data['error_norms']))
            if max_err > 5e-5:
                print("max error of", str(max_err), "is too large")
                return 1


    return 0


if __name__ == "__main__":
    sys.exit(main())
