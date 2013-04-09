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
    for f in glob.glob("results/*"):
        os.remove(f)

    # Build + run
    subp.check_call(['make', 'check', '-k'])

    # Plot trace file
    trace = sp.loadtxt("results/trace.dat", skiprows=1)
    # VARIABLES="time","u<SUB>FE</SUB>","u<SUB>exact</SUB>","norm of error","norm of solution"

    plt.plot(trace[:,0], trace[:,2], '-', label="Exact")
    plt.plot(trace[:,0], trace[:,1], 'kx', label="Approx")
    plt.plot(trace[:,0], trace[:,5], 'ro', label="dt")
    plt.legend()
    plt.show()


    return 0

# If this script is run from a shell then run main() and return the result.
if __name__ == "__main__":
    sys.exit(main())
