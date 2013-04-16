#!/usr/bin/env python

import subprocess as subp
from multiprocessing import Pool
import sys
import argparse
import os
import os.path
import multiprocessing
import itertools as it
import glob

import scipy as sp
import matplotlib.pyplot as plt

from functools import partial as pt

def main():
    """

    """

    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('--no-rerun', action='store_false', dest='rerun',
                        help="Don't delete and regenerate results")
    args = parser.parse_args()

    print args

    if args.rerun:
        # Remove old files
        for f in glob.glob("result_midpoint/*"):
            os.remove(f)

        # Build + run
        subp.check_call(['make', '-k'])
        subp.check_call(['./midpoint_llg', '-tmax', '1.0',
                         '-outdir', 'result_midpoint'])

    # Plot midpoint's trace file
    trace = sp.loadtxt("result_midpoint/trace", skiprows=1)
    # plt.plot(trace[:,0], trace[:,2], '-', label="Exact")
    # plt.plot(trace[:,0], trace[:,1], 'kx', label="Approx")
    plt.plot(trace[:,0], trace[:,2], label="midpoint Error mean")
    plt.plot(trace[:,0], trace[:,1], '+', label="midpoint dt")

    # Finish off plots
    plt.semilogy()
    plt.plot(trace[:,0], [1e-3]*len(trace[:,0]))
    plt.legend()
    plt.show()


    return 0

# If this script is run from a shell then run main() and return the result.
if __name__ == "__main__":
    sys.exit(main())
