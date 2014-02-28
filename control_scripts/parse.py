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


# Plotting functions
# ============================================================

def latex_safe(string):
    """Strip out characters that will interfere with latex.
    * _ -> space
    * ??
    """
    return string.replace("_", " ")


def make_axis_label(p, op=None):
    """Create a latex safe axis label.
    """
    if op is not None:
        return latex_safe(op.__name__ + " of " + '"' + p + '"')
    else:
        return latex_safe(p)


def plot_vs_thing(xthing, data, plot_values,
                  operations_on_values=None,
                  labels=None,
                  y_axis_lims=None):
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
        dt_label = '-dt'
    else:
        dt_label = '-tol'


    if labels is None:
        labels = ['-ref', '-ts', dt_label]
    else:
        labels = ['-ref', '-ts', dt_label] + labels

    for axes, p, op in zip(axesarray, plot_values, operations_on_values):

        for d in data:

            name = " ".join([str(d[l]) for l in labels])

            initial_vals = d.get(p)
            if initial_vals is not None:
                if op is not None:
                    vals = list(map(op, d[p]))

                else:
                    vals = d[p]


                axes.plot(d[xthing], vals, label=name)
            else:
                sys.stderr.write("Not plotting " + p
                                 + " because I couldn't find the data needed.\n")

        axes.set_ylabel(make_axis_label(p, op))

    # x label only on last axis
    axesarray[-1].set_xlabel(xthing)

    # Add a single legend
    axesarray[0].legend()

    # Resize y-axes if requested
    if y_axis_lims is not None:
        for y_lim, axis in zip(y_axis_lims, axesarray):
            if y_lim is not None:
                axis.set_ylim(y_lim)

    return fig

def plot_vs_time(*args, **kwargs):
    """Plot a list of things (plot_values) against time on a single figure with
    linked time axis.
    """
    return plot_vs_thing('times', *args, **kwargs)


def plot_vs_step(*args, **kwargs):
    """Plot a list of things (plot_values) against timestep number on a
    single figure with linked time axis.
    """
    return plot_vs_thing('DocInfo_numbers', *args, **kwargs)


def my_scatter(data, x_value, y_value, x_operation=None, y_operation=None):
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

        if x_operation is not None:
            xs = [x_operation(d[x_value]) for d in fdata]
        else:
            xs = [d[x_value] for d in fdata]

        if y_operation is not None:
            ys = [y_operation(d[y_value]) for d in fdata]
        else:
            ys = [d[y_value] for d in fdata]


        # Plot
        ax.scatter(xs, ys, s=80, marker=next(symbols), c=next(colours),
                   label=str(ts))

    # Label
    ax.set_xlabel(make_axis_label(x_value, x_operation))
    ax.set_ylabel(make_axis_label(y_value, y_operation))

    # try:
    #     # log?
    ax.set_yscale('log')
    #     ax.set_xscale('log')
    # except:
    #     return None

    # Pick the right axis scales
    ax.axis('auto')

    # lims = ax.xaxis.get_data_interval()
    # xs = sp.linspace(lims[0], lims[1])
    # ax.plot(xs, map(lambda x: x, xs),'-', label="x")
    # ax.plot(xs, map(lambda x: x**2, xs),'k-', label="x^2")
    # ax.plot(xs, map(lambda x: x**3, xs),'r-', label="x^3")
    # ax.plot(xs, map(lambda x: x**4, xs),'b-', label="x^3")

    ax.legend()



    return fig


def multi_plot(data, keys_to_split_on, plot_function):
    """Split data into sub-sets (based on keys_to_split_on) and plot each
    subset as specified by (partially evaluated) function given.
    """

    # Divide into separate data sets
    split_data = mm.split_up_stuff(data, keys_to_split_on)

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

    parser.add_argument('--split', '-s', action='append',
                        help="Split into different plots for different values of these keys, for keys begining with dash specify them as: `-s='-dt'` to avoid issues with `-` being read as a new argument.")

    args = parser.parse_args()

    if args.plots is None:
        print("No plots requested, so plotting magnetisations")
        args.plots = ['m']

    if (args.dir is None) or (args.dir == []):
        print("No directories given, so parsing ./results")
        args.dir = ["results"]

    if args.split is None:
        args.split = ['mesh', 'h-app', 'initial-m', 'mag-params', 'scale']


    if args.label is not None:
        for i, l in enumerate(args.label):
            args.label[i] = "-" + l

    # Main function
    # ============================================================

    # Get the results that aren't just empty
    really_all_results = mm.parse_parameter_sweep(args.dir)
    all_results = [d for d in really_all_results if d is not None]
    print(len(all_results), "data sets out of", len(really_all_results), "used",
          "(any others didn't have enough time steps finished).")

    print("Splitting plots based on values of", args.split)

    # Print if needed
    if args.print_data:
        pprint(all_results)


    # Do actual plots
    # ============================================================


    # Plot error norm vs time
    if 'err' in args.plots:
        plot_errors = par(plot_vs_time, plot_values=['error_norms','dts', 'trace_values'],
                          labels=args.label)
        multi_plot(all_results, args.split, plot_errors)


    # Plot m averages vs time
    if 'm' in args.plots:
        plot_m_averages = par(plot_vs_time,
                              plot_values=['mean_mxs','mean_mys','mean_mzs','dts',
                                           'h_applied_first_element'],
                              labels=args.label,
                              y_axis_lims=[[-1,1], [-1,1], [-1,1], None, None])
        multi_plot(all_results, args.split, plot_m_averages)


    if 'trace' in args.plots:
        plot_traces = par(plot_vs_time,
                          plot_values=['trace_values', 'dts'],
                          labels=args.label)
        multi_plot(all_results, args.split, plot_traces)

    if 'energy' in args.plots:
        plot_traces = par(plot_vs_time,
                          plot_values=[ 'exchange_energy',
                                        'zeeman_energy',
                                        'crystalline_anisotropy_energy',
                                        'magnetostatic_energy'],
                          labels=args.label)
        multi_plot(all_results, args.split, plot_traces)

    if 'newt' in args.plots:
        plot_newton_iters = par(plot_vs_time,
                                plot_values=['n_newton_iters','dts'],
                                labels=args.label)
        multi_plot(all_results, args.split, plot_newton_iters)

    if 'soltimes' in args.plots:
        plot_sol_time_averages = par(plot_vs_time,
                                plot_values=['solver_times','jacobian_setup_times'],
                                operations_on_values=[sp.mean, sp.mean],
                                labels=args.label)
        multi_plot(all_results, args.split, plot_sol_time_averages)

    # Plot |m| error vs time
    if 'ml' in args.plots:
        plot_ml_error_vs_time = par(plot_vs_time,
                                    plot_values=['m_length_error_means', 'dts'],
                                    labels=args.label)
        multi_plot(all_results, args.split, plot_ml_error_vs_time)


    # Plot solver iterations vs steps
    if 'its' in args.plots:
        plot_iters_step = par(plot_vs_time, plot_values=['dts', 'n_solver_iters'],
                              operations_on_values=[identity, sp.mean],
                              labels=args.label)
        multi_plot(all_results, args.split, plot_iters_step)


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
        multi_plot(all_results, args.split, plot_damping_errors)


    if 'wc-time' in args.plots:
        plot_wall_time_vs_time = \
          par(plot_vs_time,
              plot_values=['unix_timestamp_diffs', 'dts'],
              labels=args.label)

        multi_plot(all_results, args.split, plot_wall_time_vs_time)


    plt.show()

    return 0


# If this script is run from a shell then run main() and return the result.
if __name__ == "__main__":
    sys.exit(main())
