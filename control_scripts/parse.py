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

import itertools as it
import functools as ft
import scipy as sp
import matplotlib.pyplot as plt

# Imports for specific functions
from functools import partial as pt
from os.path import join as pjoin
from pprint import pprint

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


def parse_trace_file(filename):
    """Read data from a trace file into a named numpy array.

    Assuming the file is in the format

        col_Title1 col_title2 ...
        col1_value col2_value ...
        another_col1_value another_col2_value ...
        .                .
        .                .
        .                .

    """
    # ??ds input a dict of cols : names (past the first few)?
    # ??ds use structured arrays?

    # Load from the file. don't read first line (titles), only load some
    # columns (the ones I want..) and transpose the array for easy
    # extraction into seperate vectors.
    trace = sp.loadtxt(filename, skiprows=1, usecols = (1,2,3,4,21,22),
                       unpack=True, ndmin=2)

    # Unpack into separate vectors
    times, dts, err_norms, newton_iters, m_length_mean, m_length_stddev = trace

    return times, dts, err_norms, newton_iters, m_length_mean, m_length_stddev


def parse_run(results_folder):

    # If info file doesn't exist then it probably hasn't run yet...
    try:
        d = parse_info_file(pjoin(results_folder, "info"))
    except IOError:
        return {}

    times, dts, err_norms, newton_iters, mle_mean, mle_stddev \
      = parse_trace_file(pjoin(results_folder, "trace"))

    d['nsteps'] = len(times)

    d['mean_dt'] = sp.mean(dts)
    d['min_dt'] = min(dts)
    d['max_dt'] = max(dts)
    d['final_dt'] = float(dts[-1])

    d['mean_err_norm'] = sp.mean(err_norms)
    d['max_err_norm'] = max(err_norms)
    d['final_err_norm'] = err_norms[-1]

    d['mean_ml_error'] = sp.mean(mle_mean)
    d['max_ml_error'] = max(mle_mean)
    d['final_ml_error'] = mle_mean[-1]

    d['mean_newton_iters'] = sp.mean(newton_iters)
    d['max_newton_iters'] = max(newton_iters)

    # If there is a "FAILED" file then something didn't work
    d['failed'] = os.path.isfile(pjoin(results_folder, 'FAILED'))

    return d


def parse_parameter_sweep(root_dir):
    """Get a list of dictionaries of results in directories (recursively)
    contained in root_dir.
    """

    results = []
    for root, dirs, filenames in os.walk(root_dir):
        for d in dirs:
            results.append(parse_run(pjoin(root, d)))

    return results


def existing_items_match(small_dict, full_dict):
    """True if all entries that exist in the 'small dict' have the same
    value as that entry in the 'full dict'.

    Surely there should be a function for this already?
    """
    for key in small_dict:
        if full_dict[key] != small_dict[key]:
            return False

    return True


def split_up_stuff(big_dict_list):

    keys_to_split = ['initial_m', 'mesh', 'h_app', 'time_stepper', 'tmax']
    parameter_sets = [{k: bigdict[k] for k in keys_to_split}
                      for bigdict in big_dict_list]

    # use a set to get a unique list of parameter sets. Dictionaries cannot
    # be put into sets so we have to go via a tuple, i.e. list(dicts) ->
    # list(tuples(tuples)) -> set(tuples(tuples)) -> list(dicts).
    unique_parameter_sets = map(dict, set([tuple(sorted(d.items())) for d in parameter_sets]))

    newlist = []
    for test_dict in unique_parameter_sets:
        newlist.append([d for d in big_dict_list if existing_items_match(test_dict, d)])

    return newlist


def next_square_number(start):

    print("**", start)

    k = int(sp.sqrt(start))
    while k**2 < start:
        k += 1
        print(k)

    return k


def plot(data, quantity_name, use_tol=True):
    """P

    """
    if use_tol:
        y_axis_data = 'tol'
    else:
        y_axis_data = 'initial_dt'

    p = split_up_stuff(data)

    # Find out the range of the 'refinement' values
    refines = [d['refinement'] for sublist in p for d in sublist]

    # Decide if we should plot logs or normal by comparing the max and min
    # values.
    values = [d[quantity_name] for sublist in p for d in sublist]
    if max(values) > 100 * min(values):
        normaliser = sp.log10
        normaliser.label = "$\log_{10}$ of "
    else:
        # identity function
        normaliser.label = ""
        normaliser = lambda x:x

    # Make an array of subplots to put our data into
    N = len(p)
    subplt_size = next_square_number(N)
    fig, axarr = plt.subplots\
      (subplt_size, subplt_size,
       sharey=True, sharex=True, # share labels
       subplot_kw={'xlim' : (min(refines)+0.5, max(refines)+0.5),
                   'yscale' : 'log'})
    axarr = axarr.flatten()

    # Plot the data
    for ax, data_set in zip(axarr, p):
        refs = [d['refinement'] for d in data_set]
        dts = [d[y_axis_data] for d in data_set]
        vals = [normaliser(d[quantity_name]) for d in data_set]

        im = ax.scatter(refs, dts, c=vals, s=80,
                        vmin=normaliser(min(values)),
                        vmax=normaliser(max(values)))

        # Only put ticks at integer points
        ax.get_xaxis().set_major_locator(plt.MaxNLocator(integer=True))

        # Try to fit all the (other) parameters in the title
        d = data_set[0]
        ax.set_title(str(d['initial_m']) + " " + str(d['mesh']) + "\n"
                     + str(d['h_app']) + " " + str(d['time_stepper']),
                     fontsize=10)

    # Blank out any spare spaces we have left over
    for ax in axarr[N:]:
        ax.axis('off')

    fig.suptitle(quantity_name + ' for each data set')

    # "Fix" the layout
    plt.tight_layout()

    # Add a colorbar (in its own axis to make it big)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(im, cax=cbar_ax) #??ds only accurate for last subplot atm

    return fig

def main():
    """

    """


    # Parse inputs
    # ============================================================

    parser = argparse.ArgumentParser(description=main.__doc__,

    # Don't mess up my formating in the help message
    formatter_class=argparse.RawDescriptionHelpFormatter)

    args = parser.parse_args()


    # Main function
    # ============================================================


    r_all = parse_parameter_sweep("../experiments/midpoint_sweeps_2/")

    # pprint(r_all)

    r = filter(lambda d: (not d['failed']) and d['initial_m'] == 'smoothly_varying', r_all)

    pprint(len(r_all))

    dt, refine, err, mle, nsteps = ([d['mean_dt'] for d in r],
                                    [int(d['refinement']) for d in r],
                                    [d['mean_err_norm'] for d in r],
                                    [d['mean_ml_error'] for d in r],
                                    [d['nsteps'] for d in r])

    # print(dt)
    # print(refine)
    print(mle)

    # Partition data on timestepper, initial m, mesh, field

    # plot scatter of "label" vs refines with [quantity] as size (or log of
    # quantity). Color as well?

    # Replace failed runs by crosses instead of circles

    # plt.scatter(dt, refine, c=mle, s=[sp.log(m*1e6) * 50 for m in mle])
    # plt.colorbar()

    # plot(r_all, 'mean_ml_error')
    plot(r_all, 'nsteps')
    # plot(r_all, 'mean_newton_iters')
    plt.show()

    return 0


# If this script is run from a shell then run main() and return the result.
if __name__ == "__main__":
    sys.exit(main())
