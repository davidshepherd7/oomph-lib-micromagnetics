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
from functools import partial as par
from os.path import join as pjoin
from pprint import pprint


# Parsing functions
# ============================================================

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
    # extraction into seperate vectors. If there aren't enough columns
    # ignore the magnetics stuff (it's probably unsteady heat)
    try:
        return sp.loadtxt(filename, skiprows=1, usecols = (1,2,3,4,13,14,21,22,23,24,25,26),
                          unpack=True, ndmin=2)

    except IndexError:
        return sp.loadtxt(filename, skiprows=1, usecols = (1,2,3,4,13,14),
                          unpack=True, ndmin=2)


def parse_run(results_folder):

    # If info file doesn't exist then it probably hasn't run yet...
    try:
        d = parse_info_file(pjoin(results_folder, "info"))
    except IOError:
        return None

    temp = parse_trace_file(pjoin(results_folder, "trace"))

    # Unpack, only get magnetics stuff if it's there...
    try:
        (times, dts, err_norms, newton_iters, lte_est_norms, trace_values,
         mle_mean, mle_stddev,
         max_angle_variation, mxs, mys, mzs) = temp
        have_magnetics_data = True
    except ValueError:
        (times, dts, err_norms, newton_iters, lte_est_norms, trace_values) = temp
        print("Not enough fields in trace to be a magnetics driver, assuming it's not")
        have_magnetics_data = False

    # If there's only one time point then this run failed immediately and
    # we can't calculate anything interesting.
    if len(times) == 1:
        return None

    d['times'] = times
    d['dts'] = dts
    d['error_norms'] = err_norms

    d['nsteps'] = len(times[1:])

    d['mean_dt'] = sp.mean(dts[1:])
    d['min_dt'] = min(dts[1:])
    d['max_dt'] = max(dts[1:])
    d['final_dt'] = float(dts[-1])

    d['mean_err_norm'] = sp.mean(err_norms[1:])
    d['max_err_norm'] = max(err_norms[1:])
    d['final_err_norm'] = err_norms[-1]

    d['mean_newton_iters'] = sp.mean(newton_iters[1:])
    d['max_newton_iters'] = max(newton_iters[1:])

    d['lte_est_norms'] = lte_est_norms
    d['trace_values'] = trace_values

    # If there is a "FAILED" file then something didn't work
    d['failed'] = os.path.isfile(pjoin(results_folder, 'FAILED'))


    if have_magnetics_data:
        d['m_length_error_means'] = mle_mean
        d['max_angle_variation'] = max_angle_variation

        # Magnetisations
        d['mxs'] = mxs
        d['mys'] = mys
        d['mzs'] = mzs

    return d


def parse_parameter_sweep(root_dir):
    """Get a list of dictionaries of results in directories (recursively)
    contained in root_dir.
    """

    # Recursively parse directories inside the root_dir
    results = []
    for root, dirs, filenames in os.walk(root_dir):
        for d in dirs:
            print("Parsing", pjoin(root, d))
            results.append(parse_run(pjoin(root, d)))

    # Also parse the root directory if it contains results (i.e. at least
    # an info file)
    if os.path.isfile(pjoin(root, "info")):
        print("Parsing", root)
        results.append(parse_run(root))

    return results


# Other helper functions
# ============================================================


def next_square_number(start):
    k = int(sp.sqrt(start))
    while k**2 < start:
        k += 1
    return k


def existing_items_match(small_dict, full_dict):
    """True if all entries that exist in the 'small dict' have the same
    value as that entry in the 'full dict'.

    Surely there should be a function for this already?
    """
    for key in small_dict:
        if full_dict[key] != small_dict[key]:
            return False

    return True


def split_up_stuff(big_dict_list, keys_to_split_on=None):

    if keys_to_split_on is None:
       keys_to_split_on = ['initial_m', 'mesh', 'h_app', 'time_stepper', 'tmax']

    parameter_sets = [{k: bigdict[k] for k in keys_to_split_on}
                      for bigdict in big_dict_list]

    # Use a set to get a unique list of parameter sets. Dictionaries cannot
    # be put into sets so we have to go via a tuple, i.e. list(dicts) ->
    # list(tuples(tuples)) -> set(tuples(tuples)) -> unique list(dicts).
    unique_parameter_sets = map(dict, set([tuple(sorted(d.items())) for d in parameter_sets]))

    newlist = []
    for test_dict in unique_parameter_sets:
        newlist.append([d for d in big_dict_list if existing_items_match(test_dict, d)])

    return newlist


# Plotting functions
# ============================================================


def multi_scatter_plots(data, quantity_name, y_axis_data='tol'):
    """P

    """

    if 'err_norm' in quantity_name:
        # Exclude cases with no error norm given
        data = [d for d in data if d['max_err_norm'] != -1]

    # Chop the data into a list of lists (of dicts). Each list contains
    # data for the same parameters with different tol/refines.
    p = split_up_stuff(data)

    # Decide if we should plot on a log scale by comparing the max and min
    # values.
    values = [d[quantity_name] for sublist in p for d in sublist]
    print(quantity_name)
    if max(values) > 100 * min(values):
        normaliser = sp.log10
        normaliser.label = "$\log_{10}$ of "
    else:
        # identity function
        normaliser = lambda x:x
        normaliser.label = ""

    # Make an array of subplots to put our data into
    subplt_size = next_square_number(len(p))
    refines = [d['refinement'] for sublist in p for d in sublist]
    fig, axarr = plt.subplots\
      (subplt_size, subplt_size,
       sharey=True, sharex=True, # share labels
       subplot_kw={'xlim' : (min(refines)-0.5, max(refines)+0.5),
                   'yscale' : 'log'})

    # Plot the data
    for ax, data_set in zip(axarr.flatten(), p):

        # First plot ones that worked
        refs = [d['refinement'] for d in data_set if not d['failed']]
        dts = [d[y_axis_data] for d in data_set if not d['failed']]
        vals = [normaliser(d[quantity_name]) for d in data_set if not d['failed']]
        im = ax.scatter(refs, dts, c=vals, s=80, marker = 'o',
                        vmin=normaliser(min(values)),
                        vmax=normaliser(max(values)))

        # Now plot ones that failed
        refs = [d['refinement'] for d in data_set if d['failed']]
        dts = [d[y_axis_data] for d in data_set if d['failed']]
        vals = [normaliser(d[quantity_name]) for d in data_set if d['failed']]
        ax.scatter(refs, dts, c=vals, s=80, marker = 'x',
                   vmin=normaliser(min(values)),
                   vmax=normaliser(max(values)))

        # Only put ticks at integer points on x axis
        ax.get_xaxis().set_major_locator(plt.MaxNLocator(integer=True))

        # Write the (other) parameters in the title
        d = data_set[0]
        ax.set_title(str(d['initial_m']) + " " + str(d['mesh']) + "\n"
                     + str(d['h_app']) + " " + str(d['time_stepper']),
                     fontsize=10)

    # Blank out any spare spaces we have left over
    for ax in axarr.flatten()[len(p):]:
        ax.axis('off')

    fig.suptitle(quantity_name + ' for each data set')

    # "Fix" the layout
    plt.tight_layout()

    # Add a colorbar (in its own axis to make it big)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(im, cax=cbar_ax) #??ds only accurate for last subplot atm

    return fig


def plot_vs_time(data, plot_values):
    """Plot a list of things (plot_values) against time on a single figure
    with linked time axis.
    """

    fig, axesarray = plt.subplots(len(plot_values), 1, sharex = True)

    for axes, p in zip(axesarray, plot_values):

        for d in data:
            name = str(d['tol']) + " " + str(d['refinement']) + " " + str(d['time_stepper'])
            axes.plot(d['times'], d[p], label=name)
            axes.set_ylabel(p)

        axesarray[-1].set_xlabel('time')

    # Add a single legend
    axesarray[0].legend()

    return fig


def multi_plot(data, keys_to_split_on, plot_function):
    """Split data into sub-sets (based on keys_to_split_on) and plot each
    subset as specified by (partially evaluated) function given.
    """

    # Divide into separate data sets
    split_data = split_up_stuff(data, keys_to_split_on)

    for dataset in split_data:
        # Plot this one
        fig = plot_function(dataset)

        # Make a title based on the keys which specify this data set (all
        # data in dataset have the same value for the keys in
        # keys_to_split_on so we just get it from the first one).
        fig.suptitle(' '.join([dataset[0][k] for k in keys_to_split_on]))

    return


def nsteps_vs_tol(data):

    p = split_up_stuff(data, ['initial_m', 'h_app', 'mesh'])

    for data_set in p:
        fig, axarr = plt.subplots(2, 1, sharex = True)

        fig.suptitle(data_set[0]['initial_m']+ ' ' +data_set[0]['mesh'] + ' ' +
                     data_set[0]['h_app']+ ' ' +data_set[0]['time_stepper'])

        for d in [d for d in data if d['refine'] == 2]:
            axarr[0].scatter(d['tol'], sp.mean(d['error_norms']),
                          label='tol '+ str(d['tol']) +', refine '+ str(d['refinement']))
            axarr[0].set_ylabel('error norm')
            axarr[0].legend(loc=0)

            axarr[1].scatter(d['tol'], d['nsteps'],
                          label='tol '+ str(d['tol']) +', refine '+ str(d['refinement']))
            axarr[1].set_ylabel('nsteps')
            axarr[1].set_xlabel('tol')
            axarr[1].legend(loc=0)

    return


def main():
    """

    """

    # Parse inputs
    # ============================================================

    parser = argparse.ArgumentParser(description=main.__doc__,

    # Don't mess up my formating in the help message
    formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('parsedir',
                        help='Set the directory to look for data in.')

    parser.add_argument('--print-data', '-d', action='store_true',
                        help='Pretty print data to stdout')

    args = parser.parse_args()


    # Main function
    # ============================================================

    all_results = [d for d in parse_parameter_sweep(args.parsedir) if d is not None]
    print(len(all_results), "data sets")

    if args.print_data:
        pprint(all_results)

    # y_axis_data = 'mean_dt'
    # # Plot and save pdfs of interesting quantities
    # # interesting_quantities = ['mean_ml_error', 'nsteps', 'mean_newton_iters']
    # interesting_quantities = ['m_length_error_means',
    #                           ]
    # for q in interesting_quantities:
    #     fig = multi_scatter_plots(all_results, q, y_axis_data=y_axis_data)
    #     fig.savefig(pjoin(args.rootdir, q + "_plot.pdf"),
    #                 transparent=True, bbox_inches='tight')


    keys_to_split_on = ['mesh']


    # Plot error norm vs time
    plot_errors = par(plot_vs_time, plot_values=['error_norms','dts'])
    multi_plot(all_results, keys_to_split_on, plot_errors)


    # Plots that require magnetisation data (if there is none then just
    # don't plot it).
    try:
        # Plot m averages vs time
        plot_m_averages = par(plot_vs_time, plot_values=['mxs','mys','mzs','dts'])
        multi_plot(all_results, keys_to_split_on, plot_m_averages)

        # Plot |m| error vs time
        plot_ml_error_vs_time = par(plot_vs_time, plot_values=['m_length_error_means', 'dts'])
        multi_plot(all_results, keys_to_split_on, plot_ml_error_vs_time)

    except KeyError as e:
        print("KeyError with key:", e, "ignored, not doing magnetisation plots")


    plt.show()

    return 0


# If this script is run from a shell then run main() and return the result.
if __name__ == "__main__":
    sys.exit(main())
