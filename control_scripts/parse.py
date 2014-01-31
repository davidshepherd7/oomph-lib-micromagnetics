#!/usr/bin/env python3

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


# Imports for specific functions
from functools import partial as par
from os.path import join as pjoin
from pprint import pprint

# Make sure *this* versions oomphpy is in the path (before any other
# versions in other places)
sys.path.insert(1, pjoin(os.path.dirname(__file__), "../etc"))
import oomphpy
import oomphpy.micromagnetics as mm

import itertools as it
import functools as ft
import scipy as sp
import matplotlib.pyplot as plt


# Various helper functions
# ============================================================

def identity(x):
    return x


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
    refines = [d['-ref'] for sublist in p for d in sublist]
    fig, axarr = plt.subplots\
      (subplt_size, subplt_size,
       sharey=True, sharex=True, # share labels
       subplot_kw={'xlim' : (min(refines)-0.5, max(refines)+0.5),
                   'yscale' : 'log'})

    # Plot the data
    for ax, data_set in zip(axarr.flatten(), p):

        # First plot ones that worked
        refs = [d['-ref'] for d in data_set if not d['failed']]
        dts = [d[y_axis_data] for d in data_set if not d['failed']]
        vals = [normaliser(d[quantity_name]) for d in data_set if not d['failed']]
        im = ax.scatter(refs, dts, c=vals, s=80, marker = 'o',
                        vmin=normaliser(min(values)),
                        vmax=normaliser(max(values)))

        # Now plot ones that failed
        refs = [d['-ref'] for d in data_set if d['failed']]
        dts = [d[y_axis_data] for d in data_set if d['failed']]
        vals = [normaliser(d[quantity_name]) for d in data_set if d['failed']]
        ax.scatter(refs, dts, c=vals, s=80, marker = 'x',
                   vmin=normaliser(min(values)),
                   vmax=normaliser(max(values)))

        # Only put ticks at integer points on x axis
        ax.get_xaxis().set_major_locator(plt.MaxNLocator(integer=True))

        # Write the (other) parameters in the title
        d = data_set[0]
        ax.set_title(str(d['-initial-m']) + " " + str(d['-mesh']) + "\n"
                     + str(d['-h-app']) + " " + str(d['-ts']),
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


def plot_vs_time(data, plot_values, operations_on_values=None, labels=None):
    """Plot a list of things (plot_values) against time on a single figure
    with linked time axis.
    """

    # Default operations: do nothing
    if operations_on_values is None:
        operations_on_values = [None]*len(plot_values)

    # Construct axes
    fig, axesmatrix = plt.subplots(len(plot_values), 1, sharex = True,
                                   squeeze=False)
    axesarray = list(axesmatrix.flat)

    # Is dt or tol more interesting for labels? Use dt if all tols are zero
    # (i.e. non-adaptive)
    if all([d['-tol'] == 0 for d in data]):
        dt_label = 'dt'
    else:
        dt_label = 'tol'


    if labels is None:
        labels = ['ref', 'ts', dt_label]
    else:
        labesl = ['ref', 'ts', dt_label] + labels

    for axes, p, op in zip(axesarray, plot_values, operations_on_values):

        for d in data:

            # name = str(d[dt_label]) + " " + str(d['-ref'])\
               # + " " + str(d['-ts'])\
               # + " " + str(d.get('-decoupled-ms') == "1")

            name = " ".join([str(d["-"+l]) for l in labels])

            if op is not None:
                vals = map(op, d[p])
            else:
                vals = d[p]

            axes.plot(d['times'], vals, label=name)

            if op is not None:
                axes.set_ylabel(latex_safe(op.__name__ + " of " + p))
            else:
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
            name = str(d['-tol']) + " " + str(d['-ref']) + " " \
              + str(d['-ts'])

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
    for ts in set([d['-ts'] for d in data]):

        fdata = [d for d in data if d['-ts'] == ts]

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


        labels = []
        for k in keys_to_split_on:
            try:
                this_str = str(dataset[0][k])
                labels.append(this_str)

            # Ignore keys that don't exist
            except KeyError:
                pass

        # Make a title based on the keys which specify this data set (all
        # data in dataset have the same value for the keys in
        # keys_to_split_on so we just get it from the first one).
        fig.suptitle(' '.join(labels))

    return


def nsteps_vs_tol(data):

    p = split_up_stuff(data, ['-initial-m', '-h-app', '-mesh'])

    for data_set in p:
        fig, axarr = plt.subplots(2, 1, sharex = True)

        fig.suptitle(data_set[0]['-initial-m']+ ' ' +data_set[0]['-mesh'] + ' ' +
                     data_set[0]['-h-app']+ ' ' +data_set[0]['-ts'])

        for d in [d for d in data if d['-ref'] == 2]:
            axarr[0].scatter(d['-tol'], sp.mean(d['error_norms']),
                          label='-tol '+ str(d['-tol']) +', refine '+ str(d['-ref']))
            axarr[0].set_ylabel('error norm')
            axarr[0].legend(loc=0)

            axarr[1].scatter(d['-tol'], d['nsteps'],
                          label='-tol '+ str(d['-tol']) +', refine '+ str(d['-ref']))
            axarr[1].set_ylabel('nsteps')
            axarr[1].set_xlabel('-tol')
            axarr[1].legend(loc=0)

    return


def iterations_vs_dt(data):

    # Split up into separate data sets for each preconditioner
    split_data = split_up_stuff(data, ['-prec'])

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
                      c = colours.next())
                      # label=data_prec[0]['preconditioner_name'].split('-')[3])

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

    parser.add_argument('--dir', '-d', action='append',
                        help='Set the directory to look for data in (default "results").')

    parser.add_argument('--print-data', action='store_true',
                        help='Pretty print data to stdout')

    parser.add_argument('--plots', '-p', action='append',
                        help='choose what to plot (default "m")')

    parser.add_argument('--label', '-l', action='append',
                        help='Add addtional labels to line')

    args = parser.parse_args()

    if args.plots is None:
        print("No plots requested, so plotting magnetisations")
        args.plots = ['m']

    if (args.dir is None) or (args.dir == []):
        print("No directories given, so parsing ./results")
        args.dir = ["results"]

    # Main function
    # ============================================================

    # Get the results that aren't just empty
    really_all_results = mm.parse_parameter_sweep(args.dir)
    all_results = [d for d in really_all_results if d is not None]
    print(len(all_results), "data sets out of", len(really_all_results), "used",
          "(any others didn't have enough time steps finished).")

    # Print if needed
    if args.print_data:
        pprint(all_results)

    # Good default keys to split..
    keys_to_split_on = ['-mesh', '-h-app', '-initial-m', '-mag-params', '-scale']


    # Do actual plots
    # ============================================================


    # Plot error norm vs time
    if 'err' in args.plots:
        plot_errors = par(plot_vs_time, plot_values=['error_norms','dts', 'trace_values'],
                          labels=args.label)
        multi_plot(all_results, keys_to_split_on, plot_errors)


    # Plot m averages vs time
    if 'm' in args.plots:
        plot_m_averages = par(plot_vs_time,
                              plot_values=['mean_mxs','mean_mys','mean_mzs','dts'],
                              labels=args.label)
        multi_plot(all_results, keys_to_split_on, plot_m_averages)


    if 'newt' in args.plots:
        plot_newton_iters = par(plot_vs_time,
                                plot_values=['n_newton_iters','dts'],
                                labels=args.label)
        multi_plot(all_results, keys_to_split_on, plot_newton_iters)

    if 'soltimes' in args.plots:
        plot_sol_time_averages = par(plot_vs_time,
                                plot_values=['solver_times','jacobian_setup_times'],
                                operations_on_values=[sp.mean, sp.mean],
                                labels=args.label)
        multi_plot(all_results, keys_to_split_on, plot_sol_time_averages)

    # Plot |m| error vs time
    if 'ml' in args.plots:
        plot_ml_error_vs_time = par(plot_vs_time,
                                    plot_values=['m_length_error_means', 'dts'],
                                    labels=args.label)
        multi_plot(all_results, keys_to_split_on, plot_ml_error_vs_time)


    # Plot solver iterations vs steps
    if 'its' in args.plots:
        multi_plot(all_results, keys_to_split_on, iterations_vs_dt)

        plot_iters_step = par(plot_vs_time, plot_values=['dts', 'n_solver_iters'],
                              operations_on_values=[identity, sp.mean],
                              labels=args.label)
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

    if 'wc-time' in args.plots:
        plot_wall_time_vs_time = \
          par(plot_vs_time,
              plot_values=['unix_timestamp_diffs', 'dts'],
              labels=args.label)

        multi_plot(all_results, keys_to_split_on, plot_wall_time_vs_time)


    plt.show()

    return 0


# If this script is run from a shell then run main() and return the result.
if __name__ == "__main__":
    sys.exit(main())
