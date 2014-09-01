# Future proofing
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import


import string
import scipy as sp
import math
import itertools as it
import multiprocessing

from scipy.linalg import norm


# File full of generic utility functions


# General
# ============================================================

def unzip(iterable_of_iterables):
    """Inverse of zip. E.g. given a list of tuples returns a tuple of
    lists.

    To understand why: think about what * does to a list and what zip then
    does with this list.

    See http://www.shocksolution.com/2011/07/python-lists-to-tuples-and-tuples-to-lists/"""
    return list(zip(*iterable_of_iterables))


def _is_iterable(item):
    # Add checks for types you don't want to mistake as iterables here

    # Name of string base class changed in python3, try both for
    # compatability:
    try:
        isstring = isinstance(item, basestring)
    except NameError:
        isstring = isinstance(item, str)

    if isstring:
        return False

    # Fake an iteration.
    try:
        for x in item:
            break
    except TypeError:
        return False

    return True


# Parallelism
# ============================================================

def parallel_map(function, *args, serial_mode=False, **kwargs):
    """Run function with args in parallel using all available cores.

    Set serial_mode=True for debugging (uses normal map)
    """
    # For debugging we often need to run in serial (to get useful stack
    # traces).
    if serial_mode:
        results = map(function, *args)

    else:
        # Run in all parameter sets in parallel
        results = multiprocessing.Pool(**kwargs).starmap(function, unzip(args), 1)

    return list(results)

# Lists
# ============================================================

def _maybe_recursive_convert_to_array(arr):
    """Convert lists to arrays, lists of lists to arrays of arrays etc. Also
    works for any iterable. Non-iterable types are left alone.
    """
    if _is_iterable(arr):
        return sp.array([_maybe_recursive_convert_to_array(a) for a in arr])
    else:
        return arr


def differences(arr):
    """Get an array containing the differences between consecutive values of an
    array. First value is None.
    """
    return [None] + [a - b for a, b in zip(arr[1:], arr[:-1])]


def split_up_stuff(big_dict_list, keys_to_split_on=None):
    """Split a list of dicts of data into multiple lists of dicts of
    data. Split into lists based on entries in the iterable
    keys_to_split_on.
    """

    if keys_to_split_on is None:
       keys_to_split_on = []

    parameter_sets = []
    for bigdict in big_dict_list:
        thisdict = {}
        for k in keys_to_split_on:

            # If key does not exist then ignore it
            try:
                thisdict[k] = bigdict[k]
            except KeyError:
                pass

        parameter_sets.append(thisdict)


    # Use a set to get a unique list of parameter sets. Dictionaries cannot
    # be put into sets so we have to go via a tuple, i.e. list(dicts) ->
    # list(tuples(tuples)) -> set(tuples(tuples)) -> unique list(dicts).
    unique_parameter_sets = map(dict, set([tuple(sorted(d.items()))
                                           for d in parameter_sets]))

    newlist = []
    for test_dict in unique_parameter_sets:
        newlist.append([d for d in big_dict_list
                        if existing_items_match(test_dict, d)])

    return newlist


# Dictionaries
# ============================================================

def merged_dict(dict1, dict2):
    """Create a new dictionary by merging two given dictionaries.
    """
    return dict(it.chain(dict1.items(), dict2.items()))


def find_varying_value_keys(parameter_argdicts):
    """Find out which keys have values that vary between dicts

    All dicts must have the same keys.
    """

    varying_keys = set()
    for argdict in parameter_argdicts[1:]:
        for k, v in argdict.items():
            if v != parameter_argdicts[0][k]:
                varying_keys.add(k)

    return list(varying_keys)


def existing_items_match(small_dict, full_dict):
    """True if all entries that exist in the 'small dict' have the same
    value as that entry in the 'full dict'.

    Surely there should be a function for this already?
    """
    for key in small_dict:
        if full_dict.get(key) != small_dict.get(key):
            return False

    return True


# Filenames
# ============================================================


def safe_filename(unsafe):
    """Convert a string to a "safe string" which contains only ascii letters,
    numbers and underscores."""
    valid_chars = "_" + string.ascii_letters + string.digits
    return ''.join([c for c in unsafe if c in valid_chars])


def latex_escape(s):
    """Escape all characters that latex will cry about.
    """
    s = s.replace(r'{', r'\{')
    s = s.replace(r'}', r'\}')
    s = s.replace(r'&', r'\&')
    s = s.replace(r'%', r'\%')
    s = s.replace(r'$', r'\$')
    s = s.replace(r'#', r'\#')
    s = s.replace(r'_', r'\_')
    s = s.replace(r'^', r'\^{}')

    # Can't handle backslashes... ?

    return s


# Testing helpers
# ============================================================

def almost_equal(a, b, tol=1e-9):
    return abs(a - b) < tol


def abs_list_diff(list_a, list_b):
    return [abs(a - b) for a, b in zip(list_a, list_b)]


def list_almost_zero(list_x, tol=1e-9):
    return max(list_x) < tol


def list_almost_equal(list_a, list_b, tol=1e-9):
    return list_almost_zero(abs_list_diff(list_a, list_b), tol)


def same_order_of_magnitude(a, b, fp_zero):
    if abs(a) < fp_zero or abs(b) < fp_zero:
        return abs(a) < fp_zero and abs(b) < fp_zero
    else:
        return (abs(sp.log10(abs(a)) - sp.log10(abs(b))) < 1)


def same_sign(a, b, fp_zero):
    """Check if two floats (or probably fine for other numbers) have the
    same sign. Throw an error on NaN values. Treat small floats as zero and
    treat zero as not having a sign.
    """
    if (a == sp.NaN) or (b == sp.NaN):
        raise ValueError("NaN(s) passed to sign comparison functions")
    elif (abs(a) < fp_zero) and (abs(b) < fp_zero):
        return True
    else:
        return math.copysign(1, a) == math.copysign(1, b)


# Some useful asserts. We explicitly use the assert command in each
# (instead of defining the asserts in terms of the bool-returning functions
# above) to get useful output from nose -d.
def assert_almost_equal(a, b, tol=1e-9):
    assert(norm(a - b) < tol)


def assert_almost_zero(a, tol=1e-9):
    assert(norm(a) < tol)


def assert_list_almost_equal(list_a, list_b, tol=1e-9):
    assert(len(list(list_a)) == len(list(list_b)))
    for a, b in zip(list_a, list_b):
        assert(norm(a - b) < tol)


def assert_list_almost_zero(values, tol=1e-9):
    for x in values:
        assert norm(x) < tol


def assert_same_sign(a, b, fp_zero=1e-9):
    if (a == sp.NaN) or (b == sp.NaN):
        raise ValueError("NaN(s) passed to sign comparison functions")
    elif (abs(a) < fp_zero) and (abs(b) < fp_zero):
        assert True
    else:
        assert math.copysign(1, a) == math.copysign(1, b)


def assert_same_order_of_magnitude(a, b, fp_zero=1e-14):
    """Check that log10(abs(.)) are nearby for a and b. If a or b is below
    fp_zero then the other is checked in the same way for closeness to
    fp_zero (after checking that it is not also below fp_zero, for safety
    with log10).
    """
    if abs(a) < fp_zero:
        assert abs(b) < fp_zero or (
            sp.log10(abs(b)) - sp.log10(abs(fp_zero)) < 1)
    if abs(b) < fp_zero:
        assert abs(a) < fp_zero or (
            sp.log10(abs(a)) - sp.log10(abs(fp_zero)) < 1)
    else:
        assert (abs(sp.log10(abs(a)) - sp.log10(abs(b))) < 1)


# Spherical polar coordinates asserts
def assert_azi_in_range(sph):
    assert(sph.azi > 0 or almost_equal(sph.azi, 0.0))
    assert(sph.azi < 2*pi or almost_equal(sph.azi, 2*pi))


def assert_polar_in_range(sph):
    assert(sph.pol >= 0 and sph.pol <= pi)
