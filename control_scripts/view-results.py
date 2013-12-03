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
    """Convert files from .dat (tecplot) to .vtu (paraview
    compatible) in parallel. Then run paraview.

    Also do some smart things with:
    * gzipped dat files
    * Incomplete soln.pvd files"""


    # Parse inputs
    # ============================================================

    parser = argparse.ArgumentParser(description=main.__doc__,

    # Don't mess up my formating in the help message
    formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-d', '--dir', default="results/",
                        help='Pick results directory, default is "results".')

    parser.add_argument('-j', default=8, dest='ncores',
                        help='Set number of cores to use for conversions.')

    parser.add_argument('--just-unzip', '-u', action="store_true",
                        help='Just unzip the dat files zipped by -z then exit.')

    parser.add_argument('--zip-dat-files', '-z', action="store_true",
                        help='Zip up the dat files after conversion into "DIR/solndatfiles.tar.bz2".')

    args = parser.parse_args()


    # Main function
    # ============================================================

    if args.just_unzip:
        print("Just unzipping the file", pjoin(args.dir, "solndatfiles.tar.bz2"))
        subp.check_call('pbzip2 -dc ' + pjoin(args.dir, "solndatfiles.tar.bz2")
                         + ' | tar x', shell=True)
        return 0

    if len(glob.glob(pjoin(args.dir, "*.vtu"))) > 0:
        print("Already converted files, just running paraview.")

    else:

        # Unzip any individually zipped files (from parse)
        zippeddatfiles = glob.glob(pjoin(args.dir, "*.dat.gz"))
        if len(zippeddatfiles) > 0:
            subp.check_call(["gunzip"] + zippeddatfiles)
            print("gunzipped the files for you")

        # Get the list of files
        datafiles = glob.glob(pjoin(args.dir, "*.dat"))
        if len(datafiles) == 0:
            error("No .dat files found in", args.dir)

        # Convert files in parallel (map is blocking, preserves order of args).
        print("Converting", len(datafiles), "files to .vtu")
        Pool().map(convert_to_vtu, datafiles)

    # zip dat files up if requested
    if args.zip_dat_files:
        datafiles = glob.glob(pjoin(args.dir, "*.dat"))
        print("Zipping up", len(datafiles), ".dat files.")
        subp.Popen("tar -c --remove-files " + ' '.join(datafiles) + " | pbzip2 -c",
                   stdout= open(pjoin(args.dir, "solndatfiles.tar.bz2"), 'w'),
                   shell=True)

        # Notes: Popen is non-blocking (so this will run in the background
        # while paraview runs), using shell=True is evil but I can't be
        # bothered figuring out python pipes...

    # Check that soln.pvd has the final two lines, if not then add them.
    if pvd_missing_closing_lines(pjoin(args.dir, 'soln.pvd')):
        with open(pjoin(args.dir, 'soln.pvd'), 'a') as pvdfile:
            pvdfile.write('</Collection>')
            pvdfile.write('</VTKFile>')

    # Run paraview
    print("Done, launching paraview.")
    subp.call(['paraview', pjoin(args.dir, 'soln.pvd')])

    # Done!
    return 0


# If this script is run from a shell then run main() and return the result.
if __name__ == "__main__":
    sys.exit(main())
