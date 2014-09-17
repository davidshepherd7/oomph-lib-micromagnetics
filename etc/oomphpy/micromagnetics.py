# Future proofing
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

# Imports from main libraries
import subprocess as subp
import os
import os.path
import ast
import sys
import itertools as it
import functools as ft
import scipy as sp
import scipy.optimize

# Imports for specific functions
from functools import partial as par
from os.path import join as pjoin
from pprint import pprint


import oomphpy.utils as utils


# Globally accessible information functions
# ============================================================

def rootdir():
    """The location of the micromagnetics root directory."""
    return os.path.abspath(pjoin(os.path.dirname(__file__), "../../"))


def driver_path():
    """The location of the main driver binary."""
    return os.path.abspath(pjoin(rootdir(), "control_scripts/driver/driver"))


# Terminal codes for colours
def green_colour(): return '\033[01;32m'
def red_colour(): return '\033[01;31m'
def end_colour(): return '\033[0m'


def exchange_length_ms(A, Ms):
    mu0 = 4 * sp.pi * 1e-7
    return sp.sqrt((2*A) / (mu0 * (Ms**2)))


# Helper functions
# ============================================================

def badprint(*args, **kwargs):
    """Print in red"""
    args = [red_colour()] + list(args) + [end_colour()]
    print(*args, **kwargs)


def okprint(*args, **kwargs):
    """Print in green"""
    args = [green_colour()] + list(args) + [end_colour()]
    print(*args, **kwargs)


def is_arg_key(key):
    """keys are from arguments if they start with a -
    """
    return key.find("-") == 0


def convergence_rate(dataset, error_norm_norm=max):
    """For a set of convergence tests calculate the rate of covergence.
    Parameter "error_norm_norm" is used to convert the list of error
    norms at each time step to a single error norm.
    """

    # Check that we have a set of convergence tests
    assert all([d['-convergence-test'] == "1" for d in dataset])

    # Check that all parameters other than -ref match
    def arguments_match(data1, data2, exclusions=['-ref', '-outdir']):
        for k, v in data1.items():
            if is_arg_key(k) and not k in exclusions:
                if v != data2[k]:
                    return False

        return True
    assert all([arguments_match(dataset[1], d) for d in dataset[1:]])


    # Extract data to be fitted
    dts, errs = utils.unzip([( sp.mean(d['dts']), error_norm_norm(d['error_norms']) )
                             for d in dataset])

    # fit the data to a 2nd order polynomial
    def f(x, a, b, c):
        return a*(x**b) + c
    parameters, covar = sp.optimize.curve_fit(f, dts, errs)
    rate = parameters[1]     # rate of convergence = power of x (second
                             # parameter).

    return rate


def split_to_comparable_groups(dicts, key_to_compare_over, key_filter=None):
    """Given a list of dicts, split into groups whose arguments (as determined
    by is_arg_key(), outdir is excluded) only differ by
    key_to_compare_over.

    Warning: O(N^2) where N is the number of arguments in the dicts.
    """

    # Default: interested in arguments other than -outdir which always
    # changes
    if key_filter is None:
        key_filter = lambda k: (is_arg_key(k) and k != "-outdir")

    # List of keys of interest (from the first dict), assumes that all
    # dicts have the same keys.
    argument_dict_keys = [k for k, _ in dicts[0].items()
                          if (key_filter(k) and k != key_to_compare_over)]

    # List of arg keys which are not constant accross dicts. n^2 in number
    # of entries in argument_dict_keys, should be fine.
    keys_which_differ = set([])
    for k in argument_dict_keys:
        for d in dicts:
            for d2 in dicts:

                # If value differs or if key is not in one of the dicts
                # then it differs.
                try:
                    if d[k] != d2[k]:
                        keys_which_differ.add(k)
                except KeyError:
                    keys_which_differ.add(k)

    print(keys_which_differ)

    comparable_groups = utils.split_up_stuff(dicts, keys_to_split_on=keys_which_differ)

    return comparable_groups


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
        real_data = utils._maybe_recursive_convert_to_array([parse_trace_entry(c) for c in col])

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
        data['unix_timestamp_diffs'] = utils.differences(data['unix_timestamp'])



    # Calculate the actual solver time (without Jacobian assembly).
    try:
        data['real_solver_times'] = data['solver_times'] - data['jacobian_setup_times']
    except KeyError:
        data['real_solver_times'] = None

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

    # If tol != 0 then it was a time adaptive run, tol is a float though so
    # need to use a tolerance.
    d['time_adaptive'] = (abs(d['-tol']) > 1e-12)

    def integrate_time_field(d, field):
        # Get values at middle of intervals
        values = d[field]
        midpoints = [(en + enp1)/2 for en, enp1 in zip(values[:-1], values[1:])]
        return sp.dot(midpoints, d['dts'][:-1])

    d['error_norm_integral'] = integrate_time_field(d, 'error_norms')

    return d


def run_failed(dir):
    """Try to check if the run in folder dir failed. Note that this only works if the driver was run through e.g.
    """
    return any([os.path.exists(pjoin(dir, f))
                for f in ["failed", "Failed", "FAILED"]])


def parse_parameter_sweep(root_dirs, skip_failed=False, **kwargs):
    """Recursively parse directories inside each root_dir, in parallel unless
    serial_mode=True is set.

    If skip_failed then skip anything for which run_failed(dir) returns
    false.
    """
    result_dirs = []
    for root_dir in root_dirs:
        for dirname, _, _ in os.walk(root_dir):
            has_info = os.path.isfile(pjoin(dirname, "info"))
            fail_skip = skip_failed and run_failed(dirname)

            if has_info and not fail_skip:
                print("Parsing", dirname)
                result_dirs.append(dirname)
            elif fail_skip:
                print("Skipping", dirname, "because it failed")

    results = utils.parallel_map(parse_run, result_dirs, **kwargs)

    return results


# Driver execution functions
# ============================================================

def argdict2list(argdict):
    """Convert a dictionary of arguments into a command line list ready to be
    run.
    """

    special_keys = ['-mpi-ncores', '-binary', '-driver']

    # Convert any keyword args into correct format for command line input.
    processed_kwargs = []
    for key, value in argdict.items():
        if key not in special_keys:
            processed_kwargs.append(str(key))
            processed_kwargs.append(str(value))
        else:
            pass

    # Sort the arguments so that the ordering is deterministic
    sorted_processed_kwargs = list(utils.sort_as_pairs(processed_kwargs))

    # If mpi_ncores is in the dict then run with mpi and that many cores,
    # otherwise don't use mpi.
    mpi_command = []
    mpi_cores = argdict.get('-mpi-ncores')
    if mpi_cores is not None:
        mpi_command = ['mpirun', '-np', str(mpi_cores)]

    # Pull out binary path
    binary_path = argdict.get('-binary')

    # Construct argument list and command
    arglist = [str(argdict['-driver'])] + sorted_processed_kwargs

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
    # ============================================================

    with open(pjoin(outdir, "run_script"), 'w') as script_file:
        script = generate_run_script(mpi_command, binary, arglist, outdir)
        script_file.write(script)

    # Make executable (for easier re-running) (octal notation 7 =
    # read/write/execute, first digit is owner, second is group, third is
    # others, 5 = read/execute)
    os.chmod(pjoin(outdir, "run_script"), 0o775)


    # Run the command itself
    # ============================================================

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



    # Append git information to info file
    # ============================================================

    # Get micromag git info
    try:
        mm_git_rev = subp.check_output(['git', 'rev-parse', '--verify', 'HEAD'],
                                       cwd=rootdir()
                                       ).decode(sys.stdout.encoding).strip()
    except subp.CalledProcessError:
        mm_git_rev = None

    # Get oomph git info
    try:
        oomph_git_rev = subp.check_output(['git', 'rev-parse', '--verify', 'HEAD'],
                                          cwd=pjoin(rootdir(), "..", "..")
                                          ).decode(sys.stdout.encoding).strip()
    except subp.CalledProcessError:
        oomph_git_rev = None

    # Write to info file
    with open(pjoin(outdir, "info"), 'a') as info_file:
        info_file.write("micromagnetics_git_revision " + str(mm_git_rev)+"\n")
        info_file.write("oomph_lib_git_revision " + str(oomph_git_rev)+"\n")


    return err_code


def generate_run_script(mpi_command, binary, arglist, outdir):

    # Strip directory
    dirname, actual_binary = os.path.split(binary)

    # Strip -outdir and it's arg
    for i, a in enumerate(arglist):
        if a == "-outdir":
            i_outdir = i
    arglist = arglist[:i_outdir] + arglist[i_outdir+2:]

    command_string = ' '.join([binary] + arglist)

    string = "#! /bin/sh\n"
    string += "# The command run.\n"
    string += "command=\""+command_string+"\"\n"
    string += "#It should have been run in the directory cd'd into below\n"
    string += "cd " + dirname + "\n"
    string += "# Output dir args have been stripped (so that things aren't overwritten).\n"
    string += "#Output went into " + outdir + ".\n"
    string += "#stdout and stderr went into \"stdout\" in that folder\n"
    string += "\n"
    string += "# Run with gdb if any args given, otherwise just run\n"
    string += "echo \"Running\" $command\n"
    string += "if [ $# -ne 0 ]; then\n"
    string += ' '.join(mpi_command) + " gdb -ex \"run\" --args $command\n"
    string += "else\n"
    string += ' '.join(mpi_command) + " $command\n"
    string += "fi\n"

    return string


def success_message(arglist, outdir):
    okprint(' '.join(arglist))


def failure_message(arglist, outdir):
    # Print failure message to screen
    badprint(pjoin(outdir, "stdout error"))
    badprint(' '.join(arglist), "FAILED")


def generate_argdicts(main_args_dict, extra_args_dicts=None):

    if extra_args_dicts is None:
        extra_args_dicts = [{}]

    # For each extra dict create the merged dict
    arg_dicts = [utils.merged_dict(main_args_dict, extra) for extra in extra_args_dicts]

    # Make a list of all combinations of arguments in each of the dicts
    return sum([product_of_argdict(a) for a in arg_dicts], [])


def product_of_argdict(args_dict):
    """Generate a complete set of combinations of parameters."""

    # Convert any single entries into lists, be smart about strings
    # hopefully...
    for k, v in args_dict.items():
        if not utils._is_iterable(v):
            args_dict[k] = [v]

    # Create the list of dictionaries
    return [dict(zip(args_dict.keys(), x))
            for x in it.product(*args_dict.values())]


# Final function for the sweep
def _run(argdict, base_outdir, varying_args, quiet=False):

    if argdict.get('-outdir') is not None:
        error("Don't specify an outdir, it will be automatically generated")

    # Construct an informative outdir name
    varying_arg_names = [str(argdict[k]) for k in sorted(varying_args)]
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


def run_sweep(args_dict, base_outdir, extra_argsets=None, **kwargs):
    """Run a parameter sweep


    args_dict is a dictionary of argument lists.

    extra_argsets is a list of argdicts for arguments which have to come
    in groups.
    """

    # Generate list of parameter sets
    if extra_argsets is not None:
        parameter_dicts = sum([generate_argdicts(args_dict, a) for a in extra_argsets], [])
    else:
        parameter_dicts = generate_argdicts(args_dict, None)

    # Find out which arguments have different values in different runs
    varying_args = utils.find_varying_value_keys(parameter_dicts)

    # Run on all args.
    out = utils.parallel_map(_run, parameter_dicts, it.repeat(base_outdir),
                             it.repeat(varying_args), **kwargs)

    # Extract err_codes etc into separate lists and force execution (just
    # in case it's still an iterator)
    err_codes, outdirs = utils.unzip(list(out))

    return err_codes, outdirs
