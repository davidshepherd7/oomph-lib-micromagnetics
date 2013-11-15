#!/usr/bin/env python

# Future proofing
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
import sympy

import operator as op
import itertools as it
import functools as ft
import scipy as sp
import matplotlib.pyplot as plt

# Imports for specific functions
from functools import partial as par
from os.path import join as pjoin
from pprint import pprint

from matplotlib.pyplot import show as pltshow
from matplotlib.pyplot import subplots

def ts2dts(ts):
    return list(it.imap(op.sub, ts[1:], ts))


def parse_info_file(filename):
    """Read data from the info file into a dictionary.

    Assuming info file is in the format

        thing_name thing_value
    """

    info_dict = {}
    with open(filename, 'r') as f:
        for line in f:
            (key, value) = line.split()
            info_dict[key] = value

    # Make sure the numbers are stored as numbers
    info_dict['initial_dt'] = float(info_dict['initial_dt'])
    info_dict['tol'] = float(info_dict['tol'])
    info_dict['refinement'] = int(info_dict['refinement'])

    return info_dict


def main():

    ts, ys, lte_est = parse_oomph_data("results/trace")

    info_dict = parse_info_file("results/info")
    exact = sympy.lambdify(sympy.var('t'), sympy.sympify(info_dict['exact_name']+ "(t)"))
    exacts = map(exact, ts)

    fig, axes = subplots(3, 1)

    # values
    axes[0].plot(ts, ys, label="est")
    axes[0].plot(ts, exacts, label="exact")
    axes[0].legend()

    # dts
    axes[1].plot(ts[1:], ts2dts(ts))


    # ltes
    axes[2].plot(ts[5:], lte_est[5:])


    # def lte(t, dt):
    #     return abs((1/6) * dt**3 * sp.cos(t))

    # # actually there're factors of dt^3 and a constant in here too
    # lte_analyticals = map(lte, ts[1:], ts2dts(ts))

    # plt.plot(ts, lte_est, label='oomph lte')
    # plt.plot(ts[1:], lte_analyticals, label='exact lte')

    pltshow()


def parse_oomph_data(filename):

    with open(filename) as f:
        lines = f.read().splitlines()

    # first line is headers
    lines = lines[1:]

    ts = [float(l.split(";")[1]) for l in lines]
    ys = [float(l.split(";")[10]) for l in lines]
    lte_est = [float(l.split(";")[9]) for l in lines]

    return ts, ys, lte_est


# If this script is run from a shell then run main() and return the result.
if __name__ == "__main__":
    sys.exit(main())
