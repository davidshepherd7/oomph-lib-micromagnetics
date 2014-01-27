#!/usr/bin/env python3

# Python 2/3 compatability
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import subprocess as subp
import sys
import argparse
import os
import os.path
import warnings

import itertools as it
import functools as ft
import scipy as sp

# My code
import oomphpy
import oomphpy.micromagnetics as mm
import oomphpy.matrices as omat

# Imports for specific functions
from functools import partial as par
from os.path import join as pjoin

def safemax(a):
    if a != []:
        return max(a)
    else:
        return 0

def real(z): return z.real
def imag(z): return z.imag
def maxabs(a): return safemax(map(abs, a))
def maxreal(a): return safemax(map(real, a))
def maximag(a): return safemax(map(imag, a))

def main():

    results_root = pjoin(mm.rootdir(), "experiments", "parameter_sweeps", "odblock-eigs")

    for _rdir in os.listdir(results_root):
        results_dir = pjoin(results_root, _rdir)

        data = mm.parse_run(results_dir)

        if data is None:
            continue

        print("\nrefinement level", data['refinement'])
        print("disable-ms?", data['disable-ms'])

        first_m_block = 2
        nblock = 5

        # Just first Jacobian for now
        with warnings.catch_warnings():
            A = omat.import_blocks(pjoin(results_dir, "J_3_1_block"), nblock)

        # Get only LLG blocks
        A = omat.sub_matrix(A, first_m_block)

        # Get only offdiagonal blocks
        A = omat.off_diagonal_blocks(A)

        evs = []
        non_sym = []
        for block in A.flat:
            if block.nnz > 0:
                ev, evec = sp.sparse.linalg.eigs(block)
                evs.append(ev)
                TT = block - block.T

                nnz = len(filter(lambda x: x> 1e-8, TT.data))
                non_sym.append(nnz)

                # if nnz > 0:
                #     print(sp.mean(block.data))
                #     print(TT)


        # print("overall abs", map(maxabs, evs))
        # print("max real part", map(maxreal, evs))
        print("max imaginary part", max(map(maximag, evs)))

        print("max nnz in B - B^T", max(non_sym))



    return 0


if __name__ == "__main__":
    sys.exit(main())
