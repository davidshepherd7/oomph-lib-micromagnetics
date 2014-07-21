#!/usr/bin/env python3

import scipy as sp
import scipy.sparse
import scipy.sparse.linalg
import matplotlib.pyplot as plt

import oomphpy
import oomphpy.matrices as omat

from functools import partial as par
from os.path import join as pjoin

results_dir = "results"

def spy2(Z, precision=0):
    """
    SPY(Z) plots the sparsity pattern of the matrix S as an image.
    """

    from matplotlib.colors import LinearSegmentedColormap

    #binary colormap min white, max black
    cmapdata = {
    'red' : ((0., 1., 1.), (1., 0., 0.)),
    'green': ((0., 1., 1.), (1., 0., 0.)),
    'blue' : ((0., 1., 1.), (1., 0., 0.))
    }
    binary = LinearSegmentedColormap('binary', cmapdata, 2)

    Z = sp.where(abs(Z.todense()) > precision, 1.0, 0.0)

    fig, axes = plt.subplots(1, 1)
    axes.imshow(Z.T, interpolation='nearest', cmap=binary)

    return fig


def print_block_shapes(blocks):
    for i in range(0, blocks.shape[0]):
        for j in range(0, blocks.shape[1]):
            print(i, j, blocks[i, j].shape)

# Load matrices
# first Jacobian of first time step
blocks = omat.import_blocks(pjoin(results_dir, "J_1_1_block"), 7)

# Get blocks of interest
F = omat.mergeblocks(blocks[2:5, 2:5])
A = blocks[0, 0].tocsr()
D = omat.mergeblocks(sp.array(blocks[0, 2:5], ndmin=2))
B = omat.mergeblocks(sp.array([[blocks[2,0]], [blocks[3,0]], [blocks[4,0]]]))

# Some info on shapes
print("F", F.shape)
print("A", A.shape)
print("B", B.shape)
print("D", D.shape)

# Calculate the schur complement
Ainv = sp.sparse.linalg.inv(A.tocsc())
Mhat = B*Ainv*D

e, _ = omat.eig(Mhat)
print(e)
print(max(abs(e)))


# # Plots
# plt.hist(sp.log10(sp.fabs(Mhat.data)))
# spy2(omat.mergeblocks(blocks), precision=1e-10)
# plt.figure()
# plt.spy(B)


plt.show()
