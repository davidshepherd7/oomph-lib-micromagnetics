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

    # Remove old files
    for f in glob.glob("results*/*"):
        os.remove(f)

    # Build + run
    subp.check_call(['make', '-k'])
    subp.check_call(['./midpoint_llg', '-tmax', '1.0',
                     '-outdir', 'result_midpoint'])

    # Plot midpoint's trace file
    trace = sp.loadtxt("results_midpoint/trace.dat", skiprows=1)
    # plt.plot(trace[:,0], trace[:,2], '-', label="Exact")
    # plt.plot(trace[:,0], trace[:,1], 'kx', label="Approx")
    plt.plot(trace[:,0], trace[:,3], label="midpoint Error norm")
    plt.plot(trace[:,0], trace[:,4], '+', label="midpoint dt")

    # Finish off plots
    plt.plot(trace[:,0], [1e-3]*len(trace[:,0]))
    plt.legend()
    plt.show()


    return 0

# If this script is run from a shell then run main() and return the result.
if __name__ == "__main__":
    sys.exit(main())
