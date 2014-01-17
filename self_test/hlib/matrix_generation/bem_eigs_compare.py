#!/usr/bin/env python

# Python 2/3 compatability
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import scipy as sp
import itertools as it
import functools as ft
import operator as op
import sys
import sympy


# Plotting
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.pyplot import subplots
from matplotlib.pyplot import show as pltshow


# and import some common functions into the global namespace
from scipy.linalg import norm
from scipy import sin, cos, tan, log, pi, sqrt, exp, mean
from math import atan2, acos
from sympy import Rational as rat
from sympy import pretty as spretty

from oomph_python import ascii2coo
from oomph_python import ascii2array


def main():

    H = sp.loadtxt("Validation/new_bem_matrix")
    O = sp.loadtxt("Validation/old_bem_matrix")

    # H = H*-1

    # eigs
    # ============================================================


    eH, _ = sp.linalg.eig(H)
    eH = sp.array(sorted(eH))

    eO, _ = sp.linalg.eig(O)
    eO = sp.array(sorted(eO))

    print(norm(eH - eO, 2))

    fig, ax = subplots(1,1)
    ax.plot(eH, label="H")
    ax.plot(eO, label="oomph")
    ax.legend()


    # values
    # ============================================================


    vH = sp.array(sorted(H.flat))
    vO = sp.array(sorted(O.flat))

    print(norm(vH - vO, 2))


    # for h, o in zip(vH, vO):
    #     print(h, o)

    fig, ax = subplots(1,1)
    ax.plot(vH, label="H")
    ax.plot(vO, label="oomph")
    ax.legend()


    # rows
    # ============================================================
    print("\n\nrows:")
    rowdiffs = []
    for (rH, rO) in zip(H, O):
        rH = sp.sort(rH)
        rO = sp.sort(rO)
        rowdiffs.append(norm(rH - rO, 2))

    print(len(filter(lambda x: x>1e-10, rowdiffs)),
          min(rowdiffs), max(rowdiffs), mean(rowdiffs))




    # cols
    # ============================================================
    print("\n\ncols:")
    coldiffs = []
    for (cH, cO) in zip(H.T, O.T):
        cH = sp.sort(cH)
        cO = sp.sort(cO)
        coldiffs.append(norm(cH - cO, 2))

    print(len(filter(lambda x: x>1e-10, coldiffs)),
          min(coldiffs), max(coldiffs), mean(coldiffs))

    # pltshow()
    return 0


if __name__ == "__main__":
    sys.exit(main())
