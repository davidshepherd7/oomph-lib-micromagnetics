# Future proofing
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

# Imports from main libraries
import subprocess as subp
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


# Globally accessible information functions
# ============================================================

def rootdir():
    """The location of the micromagnetics root directory."""
    return os.path.abspath(pjoin(os.path.dirname(__file__), "../../"))


def driver_path():
    """The location of the main driver binary."""
    return os.path.abspath(pjoin(rootdir(), "control_scripts/driver/driver"))


def boolean_flags():
    """A list of flags which are just either enabled or not.
    """
    return ['-fd-jac', "-disable-explicit-solver-optimisations"]


# Terminal codes for colours
def green_colour(): return '\033[01;32m'
def red_colour(): return '\033[01;31m'
def end_colour(): return '\033[0m'



# Helper functions
# ============================================================
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


def unzip(l):
    """Inverse of zip. E.g. given a list of tuples returns a tuple of
    lists.

    To understand why: think about what * does to a list and what zip then
    does with this list.

    See http://www.shocksolution.com/2011/07/python-lists-to-tuples-and-tuples-to-lists/"""
    return zip(*l)


def badprint(*args, **kwargs):
    """Print in red"""
    args = [red_colour()] + list(args) + [end_colour()]
    print(*args, **kwargs)


def okprint(*args, **kwargs):
    """Print in green"""
    args = [green_colour()] + list(args) + [end_colour()]
    print(*args, **kwargs)


def existing_items_match(small_dict, full_dict):
    """True if all entries that exist in the 'small dict' have the same
    value as that entry in the 'full dict'.

    Surely there should be a function for this already?
    """
    for key in small_dict:
        if full_dict.get(key) != small_dict.get(key):
            return False

    return True


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



# Parsing functions
# ============================================================
def parse_trace_entry(entry):
    """Convert trace file entries into python data.
    """
    # ast can handle most things, catch nans first though
    if entry == "nan":
        return sp.NaN
    elif entry == "-nan":
        return -sp.NaN
    else:
        try:
            return ast.literal_eval(entry)
        except ValueError as e:
            print("Failed to parse", entry)
            return sp.NaN


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

    # Didn't even manage one step
    if len(body) == 0:
        return None

    # Make an empty dict to store our data in
    data = {}

    # Read in data from body in columns and prepend the header for that
    # column.
    for header, col in zip(headers, zip(*body)):

        # Convert the data in "col" from strings to python data types
        # (e.g. floats, lists). If it's a list (or list of lists or...)
        # type then convert to an sp.array.
        real_data = _maybe_recursive_convert_to_array([parse_trace_entry(c) for c in col])

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
    if data.get('unix_timestamp') is not None:
        data['unix_timestamp_diffs'] = differences(data['unix_timestamp'])

    return data


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
    info_dict['-dt'] = float(info_dict['-dt'])
    info_dict['-tol'] = float(info_dict['-tol'])
    info_dict['-ref'] = int(info_dict['-ref'])

    return info_dict



def parse_run(results_folder):
    """Parse the info file, trace file and look for FAILED file in a folder
    and put all data into a dict with keys taken from file data.
    """

    # If info file doesn't exist then it probably hasn't run yet...
    try:
        d = parse_info_file(pjoin(results_folder, "info"))
    except IOError:
        return None

    trace_dict = parse_trace_file(pjoin(results_folder, "trace"))

    if trace_dict is None:
        return None

    # Add the data from trace file into the dict (NOTE: if any fields have
    # the same name then the trace file data will "win").
    d.update(trace_dict)

    # If there's only one time point then this run failed immediately and
    # we can't calculate anything interesting.
    if len(d['times']) == 1:
        return None

    # If there is a "FAILED" file then something didn't work
    d['failed'] = os.path.isfile(pjoin(results_folder, 'FAILED'))

    return d


def parse_parameter_sweep(root_dirs):
    """Get a list of dictionaries of results in directories (recursively)
    contained in the list of root_dirs.
    """

    results = []
    for root_dir in root_dirs:
        # Recursively parse directories inside the root_dir
        for root, dirs, filenames in os.walk(root_dir):
            for d in dirs:
                print("Parsing", pjoin(root, d))
                results.append(parse_run(pjoin(root, d)))

        # Also parse the root directory if it contains results (i.e. at least
        # an info file)
        if os.path.isfile(pjoin(root_dir, "info")):
            print("Parsing", root_dir)
            results.append(parse_run(root_dir))

    return results


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

def parallel_map(function, *args, serial_mode=False):
    """Run function with args in parallel using all available cores.
    """

    import multiprocessing

    # For debugging we often need to run in serial (to get useful stack
    # traces).
    if serial_mode:
        results = map(function, *args)

    else:
        # Run in all parameter sets in parallel
        results = multiprocessing.Pool().starmap(function, unzip(args), 1)

    return results



# Driver execution functions
# ============================================================

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


def run_driver(arglist, outdir, binary=None, mpi_command=None):
    """Run a driver, put all output in the right places, write a script giving
    the exact command used.
    """

    # Defaults: default driver and no mpi
    if binary is None:
        binary = driver_path()
    if mpi_command is None:
        mpi_command = []


    # Always use absolute path to outdir
    # ============================================================
    outdir = os.path.abspath(outdir)

    # Make sure outdir is consistent: find its arg + compare
    try:
        outdir_index = arglist.index('-outdir')
        outdir_from_arglist = os.path.abspath(arglist[outdir_index+1])
        if outdir_from_arglist != outdir:
            print("Mismatched outdirs:", outdir_from_arglist, outdir)
            return -10

        # Replace it anway to make sure it's an absolute path
        arglist[outdir_index+1] = outdir

    # If its not in there then insert it
    except ValueError:
        arglist = arglist + ['-outdir', outdir]


    # Make sure we have strings everywhere
    arglist = [str(a) for a in arglist]

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

    # Run with specified args, and in the driver folder. Put output (stdout
    # and stderr) into a file.
    with open(pjoin(outdir, "stdout"), 'w') as stdout_file:
        err_code = subp.call(mpi_command + [binary] + arglist,
                             cwd=os.path.dirname(binary),
                             stdout = stdout_file,
                             stderr = subp.STDOUT)

    # Maybe failure message into file in ouput directory
    if err_code != 0:
        with open(pjoin(outdir, "FAILED"), 'w') as fail_file:
            print('This run failed!', file=fail_file)

    return err_code


def success_message(arglist, outdir):
    print(green_colour(), ' '.join(arglist), end_colour())


def failure_message(arglist, outdir):
    # Print failure message to screen
    print(red_colour(), pjoin(outdir, "stdout:1:1: error"), end_colour())
    print(red_colour(), ' '.join(arglist), "FAILED", end_colour())



def product_of_argdict(args_dict):
    """Generate a complete set of combinations of parameters."""

    # Convert any single entries into lists, be smart about strings
    # hopefully...
    for k, v in args_dict.items():
        if not _is_iterable(v):
            args_dict[k] = [v]

    # Create the list of dictionaries
    return [dict(zip(args_dict.keys(), x))
            for x in it.product(*args_dict.values())]


# Final function for the sweep
def _run(argdict, base_outdir, varying_args, quiet=False):

    if argdict.get('-outdir') is not None:
        error("Don't specify an outdir, it will be automatically generated")

    # Construct an informative outdir name
    varying_arg_names = sorted([str(argdict[k]) for k in varying_args])
    outdir = pjoin(base_outdir, "results_" + "_".join(varying_arg_names))

    # Make directory
    cleandir(outdir)

    # run with the args given, return err code and outdir
    arglist, binary_path, mpi_command = argdict2list(argdict)
    err = run_driver(arglist, outdir, binary_path, mpi_command)

    # Do output messages
    if not quiet:
        if err != 0:
            failure_message(arglist, outdir)
        else:
            success_message(arglist, outdir)

    return err, outdir


def _run_mp(args):
    return _run(*args)


def run_sweep(args_dict, base_outdir, parallel_sweep=False):

    # Make a list of arguments that take multiple different values
    varying_args = []
    for k, v in args_dict.items():
        if _is_iterable(v) and len(v) > 1:
            varying_args.append(k)

    # Generate list of parameter sets
    parameter_dicts = product_of_argdict(args_dict)

    # Run on all args. A little hacky because Pool() can't take locally
    # defined fuctions, can't take multiple args in map and doesn't have a
    # starmap.
    if parallel_sweep:
        args = zip(parameter_dicts, it.repeat(base_outdir),
                   it.repeat(varying_args))
        out = multiprocessing.Pool().map(_run_mp, args, 1)
    else:
        out = map(_run, parameter_dicts, it.repeat(base_outdir),
                  it.repeat(varying_args))

    # Extract err_codes etc into separate lists and force execution (just
    # in case it's still an iterator)
    err_codes, outdirs = unzip(list(out))

    return err_codes, outdirs
