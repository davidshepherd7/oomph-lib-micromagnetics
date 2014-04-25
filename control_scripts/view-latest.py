#!/usr/bin/env python3


# Some python 3 compatability. With these imports most scripts should work
# in both python 2.7 and python 3.
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

# Imports from main libraries
import subprocess as subp
from multiprocessing import Pool
import multiprocessing
import sys
import argparse
import os
import os.path
import glob

import itertools as it
import functools as ft

# Imports for specific functions
from functools import partial as pt
from os.path import join as pjoin


def error(string, code=1):
    sys.stderr.write(string)
    sys.exit(code)


def convert_to_vtu(data):
    subp.check_call(['oomph-convert', str(data)],
                    stdout = open(os.devnull, 'w'))
    return


def pvd_missing_closing_lines(pvd_filename):
    """Check that the final two lines of a pvd file are </Collection> and
    </VTKFile>, as required for reading by paraview (but often not written
    if simulation crashes).
    """
    # Read the .pvd
    with open(pvd_filename, 'r') as pvdfile:
        lines = pvdfile.readlines()

    # Check the final lines
    return (lines[-2].strip(' \t\n\r') != "</Collection>"
            or lines[-1].strip(' \t\n\r') != "</VTKFile>")


def main():
    """Convert most recent soln.dat file from .dat (tecplot) to .vtu (paraview
    compatible). Then run open in paraview. """


    # Parse inputs
    # ============================================================

    parser = argparse.ArgumentParser(description=main.__doc__,

    # Don't mess up my formating in the help message
    formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-d', '--dir', default="results/",
                        help='Pick results directory, default is "results".')

    args = parser.parse_args()


    # Main function
    # ============================================================

    # Find most recent
    most_recent_soln_file = max(glob.iglob(pjoin(args.dir, '*.dat')), key=os.path.getctime)
    file_basename = os.path.splitext(os.path.basename(most_recent_soln_file))[0]

    # Convert
    convert_to_vtu(pjoin(args.dir, file_basename + ".dat"))

    # Run paraview
    print("Done, launching paraview.")
    subp.call(['paraview', os.path.abspath(pjoin(args.dir, file_basename + ".vtu"))])

    # Done!
    return 0


# If this script is run from a shell then run main() and return the result.
if __name__ == "__main__":
    sys.exit(main())
