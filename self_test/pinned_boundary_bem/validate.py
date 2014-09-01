#!/usr/bin/env python3

# Python 2/3 compatability
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import sys
import argparse
import os
import os.path
import scipy as sp

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
        "-solver" : "gmres",
        "-matrix-type" : "som",
        "-prec" : "som-main-ilu-0",
        "-mesh" : "sq_cubeoid",
        "-max-steps" : 1,
        "-h-app" : "zero",
        "-initial-m" : "xz",
        "-hlib-bem" : 0,
        "-phi1-singularity-method" : "pin_boundary",
        "-dt" : 0.001,
        '-ts' : 'bdf2',
        "-ms-method" : ["implicit", "decoupled"],
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

    # Load energies from file
    exact_energy_array = sp.loadtxt("../micromag_energy/cubeoid/cubeoid_energies", skiprows=1)
    exact_energies = {
        "exchange_energy" : exact_energy_array[0],
        "zeeman_energy" : exact_energy_array[1],
        "crystalline_anisotropy_energy" : exact_energy_array[2],
        "magnetostatic_energy" : exact_energy_array[3],
        }
    energy_keys = [k for k,_ in exact_energies.items()]

    # Check answers
    t1 = [tests.check_dicts_match(d, exact_energies, energy_keys, rtol=8, zero = 1.5e-4)
          for d in datasets]

    if ran and all(t1):
        return 0
    else:
        return 1


if __name__ == "__main__":
    sys.exit(main())
