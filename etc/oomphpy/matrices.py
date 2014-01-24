import scipy as sp
import scipy.sparse
import scipy.sparse.linalg

# For some reason numpy.object doesn't seem to be included in scipy
# (yet?). Everything else is in scipy nowadays.
import numpy

# Read in matrices/vectors from oomph formatted output
# ============================================================

def ascii2coo(filename):
    """Convert from an ascii file "filename" contatining a matrix in coordinate
    notation to scipy sparse in coordinate notation."""

    ascii_contents = scipy.loadtxt(filename,ndmin=2)
    data = ascii_contents[:,2]
    ij = [ascii_contents[:,0], ascii_contents[:,1]]
    return scipy.sparse.coo_matrix((data,ij));


def ascii2array(filename):
    """Read an oomph-lib double vector into an array.
    """
    return scipy.loadtxt(filename, usecols=[1])


def import_blocks(base_filename, nblocks):
    """Get a "matrix of sparse matrices" by reading from filenames in the format
    "base_filename_i_j" for i,j in 0:nblocks. We use coordinate list sparse
    matrices here for ease of modification. Convert to csr later."""

    block_array = scipy.zeros((nblocks,nblocks), dtype=numpy.object)
    for i in range(nblocks):
        for j in range(nblocks):
            filename = "{0}_{1}_{2}".format(base_filename,i,j)
            block_array[i,j] = ascii2coo(filename)
    return block_array;


# Functions to extract only some blocks
# ============================================================
def throw_away_blocks(blocked_matrix, condition_function):
    """If the function condition(index) is not true for the (i,j)th block then
    replace it with an empty sparse matrix of the correct size. We use
    coordinate list sparse matrices here for ease of modification. Convert to
    csr later."""

    for index, x in scipy.ndenumerate(blocked_matrix):
        if condition_function(index[0], index[1]):
            blocked_matrix[index] = scipy.sparse.coo_matrix(blocked_matrix[index].get_shape())
    return blocked_matrix


def diagonal_block_matrix(A):
    return throw_away_blocks(A, lambda i, j: i != j)


def block_upper_triangular(A):
    return throw_away_blocks(A, lambda i, j: i > j)


def lower_upper_triangular(A):
    return throw_away_blocks(A, lambda i, j: i < j)


def sub_matrix(A, first_block=0, last_block=None):
    if last_block is None:
        last_block = A.shape[0]

    def f(i, j):
        i_drop = i < first_block or i > last_block
        j_drop = j < first_block or j > last_block
        return i_drop or j_drop

    return throw_away_blocks(A, f)


def off_diagonal_blocks(A):
    return throw_away_blocks(A, lambda i, j: i ==j)


#  Other functions
# ============================================================
def mergeblocks(A):
    """To merge block matrices just use bmat(blocked_matrix) then convert to
    compressed row form ready for use in solvers."""
    return scipy.sparse.csr_matrix(scipy.sparse.bmat(A));


def print_block_sparsity(A):
    for index, x in scipy.ndenumerate(A):
        if issparse(A) and isempty(A):
            D[index] = 0
        elif issparse(A):
            D[index] = F

# ??ds when you need dense sub blocks remember to look at bsr_matrix which is perfect.
