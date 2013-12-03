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
import ast

import itertools as it
import functools as ft
import scipy as sp
import matplotlib.pyplot as plt

# Imports for specific functions
from functools import partial as par
from os.path import join as pjoin
from pprint import pprint

def identity(x):
    return x

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
        # (e.g. floats, lists).
        real_data = [ast.literal_eval(c) for c in col]

        # Put it into our dict with the header of the column as the key and
        # the rest as the data.
        data[header] = real_data


    # Remove 'dummy' columns if we had any (may have used them for padding)
    try:
        del data['dummy']
    except KeyError:
        pass

    return data


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
    if os.path.isfile(pjoin(root_dir, "info")):
        print("Parsing", root_dir)
        results.append(parse_run(root_dir))

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
    """Split a list of dicts of data into multiple lists of dicts of
    data. Split into lists based on entries in the iterable
    keys_to_split_on.
    """

    if keys_to_split_on is None:
       keys_to_split_on = []

    parameter_sets = [{k: bigdict[k] for k in keys_to_split_on}
                      for bigdict in big_dict_list]

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


def latex_safe(string):
    """Strip out characters that will interfere with latex.
    * _ -> space
    * ??
    """
    return string.replace("_", " ")


def plot_vs_time(data, plot_values):
    """Plot a list of things (plot_values) against time on a single figure
    with linked time axis.
    """

    fig, axesmatrix = plt.subplots(len(plot_values), 1, sharex = True,
                                   squeeze=False)
    axesarray = axesmatrix.flat

    for axes, p in zip(axesarray, plot_values):

        for d in data:
            name = str(d['tol']) + " " + str(d['refinement'])\
               + " " + str(d['time_stepper'])\
               + " " + str(d.get('use_implicit_ms') is not None)

            axes.plot(d['times'], d[p], label=name)
            axes.set_ylabel(latex_safe(p))

        axesarray[-1].set_xlabel('time')

    # Add a single legend
    axesarray[0].legend()

    return fig


def plot_vs_step(data, plot_values, operations=None):
    """Plot a list of things (plot_values) against timestep number on a
    single figure with linked time axis. If operations is a list of
    functions then call each function on the corresponding plot_value data
    before plotting.
    """

    if operations is None:
        operations = [identity] * len(plot_values)

    fig, axesmatrix = plt.subplots(len(plot_values), 1, sharex = True,
                                   squeeze=False)
    axesarray = axesmatrix.flat

    for axes, p, f in zip(axesarray, plot_values, operations):

        for d in data:
            name = str(d['tol']) + " " + str(d['refinement']) + " " \
              + str(d['time_stepper'])

            axes.plot([f(y) for y in d[p]], label=name)

            if f is not identity:
                axes.set_ylabel(str(f) + " of " + latex_safe(p))
            else:
                axes.set_ylabel(latex_safe(p))

        axesarray[-1].set_xlabel('time step')

    # Add a single legend
    axesarray[0].legend()

    return fig


def my_scatter(data, x_value, y_value, x_operation=sp.mean, y_operation=sp.mean):
    """Plot a scatter plot at point (x_value, y_value) for each dataset in
    data. Perform operations on data before plotting (by default take
    mean).
    """

    fig, ax = plt.subplots()

    symbols = iter(['x', 'o', '+', '^', '*'])
    colours = iter(['r', 'g', 'b', 'k', 'c'])

    # Plot each time stepper
    for ts in set([d['time_stepper'] for d in data]):

        fdata = [d for d in data if d['time_stepper'] == ts]

        xs = [x_operation(d[x_value]) for d in fdata]
        ys = [y_operation(d[y_value]) for d in fdata]

        # Plot
        ax.scatter(xs, ys, s=80, marker=symbols.next(), c=colours.next(),
                   label=str(ts))

    # Label
    ax.set_xlabel(str(x_operation) + " of " + x_value)
    ax.set_ylabel(str(y_operation) + " of " + y_value)

    try:
        # log?
        ax.set_yscale('log')
        ax.set_xscale('log')
    except:
        return None

    # Pick the right axis scales
    ax.axis('auto')

    lims = ax.xaxis.get_data_interval()
    xs = sp.linspace(lims[0], lims[1])
    ax.plot(xs, map(lambda x: x, xs),'-', label="x")
    ax.plot(xs, map(lambda x: x**2, xs),'k-', label="x^2")
    ax.plot(xs, map(lambda x: x**3, xs),'r-', label="x^3")
    ax.plot(xs, map(lambda x: x**4, xs),'b-', label="x^3")

    ax.legend()



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
        fig.suptitle(' '.join([str(dataset[0][k]) for k in keys_to_split_on]))

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


def iterations_vs_dt(data):

    # Split up into separate data sets for each preconditioner
    split_data = split_up_stuff(data, ['preconditioner_name'])

    # Create figure
    fig, axarr = plt.subplots(1, 1)

    symbols = iter(['x', 'o', '+', '^', '*'])
    colours = iter(['r', 'g', 'b', 'k', 'c'])

    # Plot as scatter vs dt
    for data_prec in split_data:

        iterations = [sp.mean([sp.mean(d2) for d2 in d['n_solver_iters'][1:]])
                      for d in data_prec]
        dts = [sp.mean(d['dts']) for d in data_prec]

        axarr.scatter(dts, iterations, marker = symbols.next(), s = 50,
                      c = colours.next(),
                      label=data_prec[0]['preconditioner_name'].split('-')[3])

        axarr.set_xlabel("dt")
        axarr.set_ylabel("N solver iterations")


    axarr.legend()

    return fig


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

    parser.add_argument('--print-data', action='store_true',
                        help='Pretty print data to stdout')

    parser.add_argument('--plots', '-p', action='append',
                        help='choose what to plot')

    args = parser.parse_args()

    if args.plots is None:
        print("No plots requested, so plotting magnetisations")
        args.plots = ['m']

    # Main function
    # ============================================================

    # Get the results that aren't just empty
    really_all_results = parse_parameter_sweep(args.parsedir)
    all_results = [d for d in really_all_results if d is not None]
    print(len(all_results), "data sets out of", len(really_all_results), "used",
          "(any others didn't have enough time steps finished).")

    # Print if needed
    if args.print_data:
        pprint(all_results)

    # Good default keys to split..
    keys_to_split_on = ['mesh', 'refinement', 'h_app', 'initial_m', 'mag_params']


    # Do actual plots
    # ============================================================


    # Plot error norm vs time
    if 'err' in args.plots:
        plot_errors = par(plot_vs_time, plot_values=['error_norms','dts', 'trace_values'])
        multi_plot(all_results, keys_to_split_on, plot_errors)


    # Plot m averages vs time
    if 'm' in args.plots:
        plot_m_averages = par(plot_vs_time,
                              plot_values=['mean_mxs','mean_mys','mean_mzs','dts'])
        multi_plot(all_results, keys_to_split_on, plot_m_averages)


    # Plot |m| error vs time
    if 'ml' in args.plots:
        plot_ml_error_vs_time = par(plot_vs_time,
                                    plot_values=['m_length_error_means', 'dts'])
        multi_plot(all_results, keys_to_split_on, plot_ml_error_vs_time)


    # Plot solver iterations vs steps
    if 'its' in args.plots:
        multi_plot(all_results, keys_to_split_on, iterations_vs_dt)

        plot_iters_step = par(plot_vs_step, plot_values=['dts', 'n_solver_iters'],
                              operations=[identity, sp.mean])
        multi_plot(all_results, keys_to_split_on, plot_iters_step)


    # Plot error in effective damping vs step size
    if 'damp' in args.plots:

        def damping_error_mean(effective_dampings, exact_damping=0.0):
            # Drop first few steps: they are wrong
            effective_dampings = effective_dampings[6:]

            # relative error
            if exact_damping != 0.0:
                e = abs(sp.array(effective_dampings) - exact_damping)/exact_damping
            else:
                e = abs(sp.array(effective_dampings))

            # take mean
            return sp.mean(e)

        plot_damping_errors = par(my_scatter, x_value='dts', y_value='effective_damping',
                                  y_operation=damping_error_mean)
        multi_plot(all_results, keys_to_split_on, plot_damping_errors)


    plt.show()

    return 0


# If this script is run from a shell then run main() and return the result.
if __name__ == "__main__":
    sys.exit(main())
