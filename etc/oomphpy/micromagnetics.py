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
import ast

import itertools as it
import functools as ft
import scipy as sp
import matplotlib.pyplot as plt

# Imports for specific functions
from functools import partial as par
from os.path import join as pjoin
from pprint import pprint

def rootdir():
    """Get the micromagnetics root directory"""
    return os.path.abspath(pjoin(os.path.dirname(__file__), "../../"))

def driver_path():
    """The location of the main driver binary"""
    return pjoin(rootdir(), "control_scripts/driver/driver")


def cleandir(dirname):
    """(Re)make a directory called dirname.
    """

    # If it exists then delete all files (won't touch subdirs though) and
    # the folder itself. This will fail if we have any subdirs (for safety).
    if os.path.isdir(dirname):
        for f in os.listdir(dirname):
            os.unlink(pjoin(dirname, f))
        os.rmdir(dirname)

    # Make the directory and any parents needed
    os.makedirs(dirname)


def _is_iterable(item):
    # Add types here you don't want to mistake as iterables
    if isinstance(item, basestring):
        return False

    # Fake an iteration.
    try:
        for x in item:
            break
    except TypeError:
        return False

    return True


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


def parse_trace_file(filename):
    """Read data from a trace file into a dict of lists keyed on the column
    header.

    Assuming the file is in the format

        col1_header col2_header ...
        col1_value1 col2_value1 ...
        col1_value2 col2_value2 ...
        .                .
        .                .
        .                .

    """
    # Convert to new deliminator with "find -name 'trace' | xargs sed -e
    # 's/ /; /g' -i" might not work for very big files due to limitations
    # in sed -i?

    # Get the file data as a list of strings
    with open(filename) as f:
        lines = f.read().splitlines()

    def mysplit(line):
        # Split on '; ' and skip the last entry (because we have a final ';
        # ' at the end of lines).
        return line.split('; ')#[:-1]

    # Split the raw strings into data fields (so we now have a list of strings)
    headers = mysplit(lines[0])
    body = [mysplit(l) for l in lines[1:]]

    # Make an empty dict to store our data in
    data = {}

    # Read in data from body in columns and prepend the header for that
    # column.
    for header, col in zip(headers, zip(*body)):

        # Convert the data in "col" from strings to python data types
        # (e.g. floats, lists). If it's a list (or list of lists or...)
        # type then convert to an sp.array.
        real_data = _maybe_recursive_convert_to_array([ast.literal_eval(c) for c in col])

        # Put it into our dict with the header of the column as the key and
        # the rest as the data.
        data[header] = real_data


    # Remove 'dummy' columns if we had any (may have used them for padding)
    try:
        del data['dummy']
    except KeyError:
        pass

    # Calculate time stamp differences (rough value for wall clock time per
    # step).
    data['unix_timestamp_diffs'] = differences(data['unix_timestamp'])

    return data


def boolean_flags():
    """A list of flags which are just either enabled or not. Inside function so
    it's global but harder to accidentally modify.
    """
    return ['-decoupled-ms', '-disable-ms', '-fd-jac']


def argdict2list(argdict):
    """Convert a dictionary of arguments into a command line list ready to be
    run.
    """

    special_keys = ['-mpi-ncores', '-binary', '-driver'] + boolean_flags()

    # Convert any keyword args into correct format for command line input.
    processed_kwargs = []
    for key, value in argdict.items():
        if key not in special_keys:
            processed_kwargs.append(str(key))
            processed_kwargs.append(str(value))

        # If it's a bool flag then either add it or don't, depending on the
        # boolean value
        elif key in boolean_flags():
            if value:
                processed_kwargs.append(str(key))
            else:
                pass


    # If mpi_ncores is in the dict then run with mpi and that many cores,
    # otherwise don't use mpi.
    mpi_command = []
    mpi_cores = argdict.get('-mpi-ncores')
    if mpi_cores is not None:
        mpi_command = ['mpirun', '-np', str(mpi_cores)]

    # Pull out binary path
    binary_path = argdict.get('-binary')

    # Construct argument list and command
    arglist = [str(argdict['-driver'])] + processed_kwargs

    return arglist, binary_path, mpi_command


def parallel_parameter_sweep(function, parameter_dictionary, serial_mode=False):
    """Run function with all combinations of parameters in parallel using
    all available cores.

    parameter_lists should be a list of lists of parameters,
    """

    import multiprocessing

    # Generate a complete set of combinations of parameters
    parameter_sets = [dict(zip(parameter_dictionary, x))
                      for x in it.product(*parameter_dictionary.values())]

    # For debugging we often need to run in serial (to get useful stack
    # traces).
    if serial_mode:
        results = map(function, parameter_sets)

    else:
        # Run in all parameter sets in parallel
        results = multiprocessing.Pool().map(function, parameter_sets, 1)

    return results


def run_driver(arglist, outdir, binary=None, mpi_command=None):
    """Run a driver, put all output in the right places, write a script giving
    the exact command used.
    """

    # Defaults: default driver and no mpi
    if binary is None:
        binary = driver_path()
    if mpi_command is None:
        mpi_command = []

    # Check that the binary exists
    if not os.path.isfile(binary):
        print("Can't find binary", binary)
        raise OSError;

    # Write the command used to a file
    with open(pjoin(outdir, "run_script"), 'w') as script_file:
        script_file.write("#!/bin/sh\n")
        script_file.write("# Command used in run\n")
        script_file.write(' '.join(mpi_command + [binary] + arglist))
        script_file.write("\n")
        script_file.write("# stdout and stderr go into "
                          + pjoin(outdir, "stdout") + "\n")

    # Run with specified args, put output (stdout and stderr) into a file.
    with open(pjoin(outdir, "stdout"), 'w') as stdout_file:
        err_code = subp.call(mpi_command + [binary] + arglist,
                             stdout = stdout_file,
                             stderr = subp.STDOUT)

    return err_code
