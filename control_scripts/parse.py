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
from shutil import copy as cp

# Make sure *this* versions oomphpy is in the path (before any other
# versions in other places)
sys.path.insert(1, pjoin(os.path.dirname(__file__), "../etc"))
import oomphpy
import oomphpy.micromagnetics as mm

import itertools as it
import functools as ft
import scipy as sp


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


def make_label(p, op=None):
    """Create a label from field name and operation.
    """
    if op is not None:
        return op.__name__ + " of " + '"' + p + '"'
    else:
        return p


def make_axis_label(*args, **kwargs):
    return latex_safe(make_label(*args, **kwargs))


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
    import matplotlib.pyplot as plt # local import needed so we can vary
                                    # the backend dynamically
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


def my_scatter(data, x_value, y_value,
               dataset_split_keys=None,
               labels=None,
               x_operation=None, y_operation=None,
               **kwargs):
    """Plot a scatter plot at point (x_value, y_value) for each dataset in
    data. Perform operations on data before plotting (by default take
    mean).
    """

     # local import needed so we can vary the backend dynamically
    import matplotlib.pyplot as plt

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

    if dataset_split_keys is None:
        dataset_split_keys = ["-damping", "-dt", "-scale"]
    else:
        labels=dataset_split_keys

    fig, ax = plt.subplots(subplot_kw=kwargs)


    a = list(it.product(['r', 'g', 'b', 'k', 'c'], ['x', 'o', '+', '^', '*']))
    symbols = it.cycle(map(''.join, a))

    split_data = mm.split_up_stuff(data, dataset_split_keys)


    # Plot each time stepper
    for fdata in split_data:

        name = " ".join([str(fdata[0][l]) for l in labels])

        if x_operation is not None:
            xs = [x_operation(d[x_value]) for d in fdata]
        else:
            xs = [d[x_value] for d in fdata]

        if y_operation is not None:
            ys = [y_operation(d[y_value]) for d in fdata]
        else:
            ys = [d[y_value] for d in fdata]

        # Plot
        ax.plot(xs, ys, next(symbols),
                markersize=10,
                label=name)

    # Label
    ax.set_xlabel(make_axis_label(x_value, x_operation))
    ax.set_ylabel(make_axis_label(y_value, y_operation))

    ax.set_xscale("log")
    ax.set_yscale("log")


    # Pick the right axis scales
    ax.axis('auto')

    ax.legend()



    return fig


def multi_plot(data, keys_to_split_on, plot_function):
    """Split data into sub-sets (based on keys_to_split_on) and plot each
    subset as specified by (partially evaluated) function given.
    """

    # Divide into separate data sets
    split_data = mm.split_up_stuff(data, keys_to_split_on)

    # Storage for figures
    figs = []

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

        figs.append(fig)

    return figs


def multi_print(data, keys_to_split_on, print_function):

    # Divide into separate data sets
    split_data = mm.split_up_stuff(data, keys_to_split_on)

    for dataset in split_data:

        # Make a title based on the keys which specify this data set (all
        # data in dataset have the same value for the keys in
        # keys_to_split_on so we just get it from the first one).
        labels = []
        for k in keys_to_split_on:
            try:
                this_str = str(dataset[0][k])
                labels.append(this_str)

            # Ignore keys that don't exist
            except KeyError:
                pass

        print(' '.join(labels))

        # Print data
        print_function(dataset)


    return


def data_print(datasets, to_print, delim="; "):

    for a, op in to_print:
        print(make_label(a, op), end=delim)

    print()

    for d in datasets:
        for a, op in to_print:
            if op is None:
                print(d[a], end=delim)
            else:
                print(op(d[a]), end=delim)

        print()

    print()


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

    parser.add_argument('--print-all-data', action='store_true',
                        help='Pretty print all data to stdout')

    parser.add_argument('--plots', '-p', action='append',
                        help='choose what to plot (default "m")')

    parser.add_argument('--print-data', '-v', action='append', default=[],
                        help='Choose values to print')

    parser.add_argument('--label', '-l', action='append',
                        help='Add additional labels to line, for keys begining with dash specify them as: `-s=-dt` to avoid issues with `-` being read as a new argument.')

    parser.add_argument('--split', '-s', action='append',
                        help="Split into different plots for different values of these keys, for keys begining with dash specify them as: `-s=-dt` to avoid issues with `-` being read as a new argument.")

    parser.add_argument('--scatter-split', '-t', action='append',
                        help="Split into different data sets in a scatter plot for different values of these keys, for keys begining with dash specify them as: `-s=-dt` to avoid issues with `-` being read as a new argument.")

    parser.add_argument('--skip-failed', action='store_true',
                        help='Skip runs which failed (dir contains file named failed)')

    parser.add_argument('--filter', '-f', action='append', default=[],
                        help="""Filter runs based on parameters. Input should be a python pair like ('-tol', 1e-3) where the first entry is the dictionary key to check and the second entry is the value it should take. Note the quotes around the key.""")

    parser.add_argument('--save-to-dir',
                        help='Save figures as pdfs into the specified folder.')

    parser.add_argument('--ssh-mode', action='store_true',
                        help="Vastly speed up operation over ssh: don't show plots and use Agg backend.")

    parser.add_argument('--save-data-to-dir', action='store_true',
                        help='Also copy trace and info files to the directory.')

    args = parser.parse_args()

    # Handle some defaults like this instead of inside argparse otherwise
    # arguments are appended to defaults rather than replacing them.
    if args.plots is None:
        print("No plots requested, so plotting magnetisations")
        args.plots = ['m']

    if (args.dir is None) or (args.dir == []):
        print("No directories given, so parsing ./results")
        args.dir = ["results"]

    if args.split is None:
        args.split = ['mesh', 'h-app', 'initial-m', 'mag-params', 'scale']

    # If we are in ssh mode then don't use a plotting type which requires
    # X11. This has to happen before we import pyplot.
    if args.ssh_mode:
        import matplotlib
        matplotlib.use('Agg')

    import matplotlib.pyplot as plt



    # Main function
    # ============================================================

    # Get the results that aren't just empty
    really_all_results = mm.parse_parameter_sweep(args.dir,
                                                  skip_failed=args.skip_failed)
    all_results = [d for d in really_all_results if d is not None]

    print(len(all_results), "data sets out of", len(really_all_results), "used",
          "(any others didn't have enough time steps finished).")

    print("Splitting plots based on values of", args.split)
    print("Labeling  plots based on values of", args.label)

    # Print if needed
    if args.print_all_data:
        pprint(all_results)

    # Filter the results based on the arguments given
    for f in args.filter:
        print("filtering with", f)
        key, value = ast.literal_eval(f)
        all_results = [d for d in all_results if d[key] == value]


    # Do actual plots
    # ============================================================

    # Storage for figures
    figs = []


    # Plot error norm vs time
    if 'err' in args.plots:
        plot_errors = par(plot_vs_time, plot_values=['error_norms','dts', 'trace_values'],
                          labels=args.label)
        newfigs = multi_plot(all_results, args.split, plot_errors)
        figs.extend(newfigs)


    # Plot m averages vs time
    if 'm' in args.plots:
        plot_m_averages = par(plot_vs_time,
                              plot_values=['mean_mxs','mean_mys','mean_mzs','dts',
                                           'h_applied_first_element'],
                              labels=args.label,
                              y_axis_lims=[[-1,1], [-1,1], [-1,1], None, None])
        newfigs = multi_plot(all_results, args.split, plot_m_averages)
        figs.extend(newfigs)


    if 'trace' in args.plots:
        plot_traces = par(plot_vs_time,
                          plot_values=['trace_values', 'dts'],
                          labels=args.label)
        newfigs = multi_plot(all_results, args.split, plot_traces)
        figs.extend(newfigs)

    if 'energy' in args.plots:
        plot_traces = par(plot_vs_time,
                          plot_values=[ 'exchange_energy',
                                        'zeeman_energy',
                                        'crystalline_anisotropy_energy',
                                        'magnetostatic_energy'],
                          labels=args.label)
        newfigs = multi_plot(all_results, args.split, plot_traces)
        figs.extend(newfigs)

    if 'newt' in args.plots:
        plot_newton_iters = par(plot_vs_time,
                                plot_values=['n_newton_iters','dts'],
                                labels=args.label)
        newfigs = multi_plot(all_results, args.split, plot_newton_iters)
        figs.extend(newfigs)

    if 'soltimes' in args.plots:
        plot_sol_time_averages = par(plot_vs_time,
                                plot_values=['solver_times','jacobian_setup_times'],
                                operations_on_values=[sp.mean, sp.mean],
                                labels=args.label)
        newfigs = multi_plot(all_results, args.split, plot_sol_time_averages)
        figs.extend(newfigs)

    # Plot |m| error vs time
    if 'ml' in args.plots:
        plot_ml_error_vs_time = par(plot_vs_time,
                                    plot_values=['m_length_error_means', 'dts'],
                                    labels=args.label)
        newfigs = multi_plot(all_results, args.split, plot_ml_error_vs_time)
        figs.extend(newfigs)

    if 'lte' in args.plots:
        plot_ml_error_vs_time = par(plot_vs_time,
                                    plot_values=['LTE_norms', 'dts'],
                                    labels=args.label)
        newfigs = multi_plot(all_results, args.split, plot_ml_error_vs_time)
        figs.extend(newfigs)


    # Plot solver iterations vs steps
    if 'its' in args.plots:
        plot_iters_step = par(plot_vs_time, plot_values=['dts', 'n_solver_iters'],
                              operations_on_values=[identity, sp.mean],
                              labels=args.label)
        newfigs = multi_plot(all_results, args.split, plot_iters_step)
        figs.extend(newfigs)


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
                                  dataset_split_keys=args.scatter_split,
                                  labels=args.label,
                                  y_operation=damping_error_mean)
        newfigs = multi_plot(all_results, args.split, plot_damping_errors)
        figs.extend(newfigs)


    if 'wc-time' in args.plots:
        plot_wall_time_vs_time = \
          par(plot_vs_time,
              plot_values=['unix_timestamp_diffs', 'dts'],
              labels=args.label)

        newfigs = multi_plot(all_results, args.split, plot_wall_time_vs_time)
        figs.extend(newfigs)


    if 'step-times' in args.plots:
        # for d in all_results:
            # print(sp.mean(d['total_step_time']))

        plot_mean_step_times_scatter = \
          par(my_scatter,
              labels=args.label,
              dataset_split_keys=args.scatter_split,
              x_value='initial_nnode',
              y_value='total_step_time',
              y_operation=sp.mean)

        newfigs = multi_plot(all_results, args.split, plot_mean_step_times_scatter)
        figs.extend(newfigs)


    if 'scatter-its' in args.plots:
        plot_mean_step_times_scatter = \
          par(my_scatter,
              labels=args.label,
              dataset_split_keys=args.scatter_split,
              x_value='initial_nnode',
              y_value='n_solver_iters',
              y_operation=lambda x:sp.mean(list(it.chain(*x))))

        newfigs = multi_plot(all_results, args.split, plot_mean_step_times_scatter)
        figs.extend(newfigs)

    if 'scatter-err-dts' in args.plots:
        plot_err_scatter = \
          par(my_scatter,
              labels=args.label,
              dataset_split_keys=args.scatter_split,
              x_value='dts',
              y_value='error_norms',
              y_operation=max,
              x_operation=sp.mean)

        newfigs = multi_plot(all_results, args.split, plot_err_scatter)
        figs.extend(newfigs)

    if 'scatter-err-tols' in args.plots:
        plot_err_scatter = \
          par(my_scatter,
              labels=args.label,
              dataset_split_keys=args.scatter_split,
              x_value='-tol',
              y_value='error_norms',
              y_operation=max)

        newfigs = multi_plot(all_results, args.split, plot_err_scatter)
        figs.extend(newfigs)

    if 'scatter-first-err' in args.plots:
        plot_err_scatter = \
          par(my_scatter,
              labels=args.label,
              dataset_split_keys=args.scatter_split,
              x_value='-dt',
              y_value='error_norms',
              y_operation=lambda y:y[1])

        newfigs = multi_plot(all_results, args.split, plot_err_scatter)
        figs.extend(newfigs)


    if 'scatter-newt' in args.plots:
        plot_mean_step_times_scatter = \
          par(my_scatter,
              labels=args.label,
              dataset_split_keys=args.scatter_split,
              x_value='initial_nnode',
              y_value='n_newton_iters',
              y_operation=sp.mean)

        newfigs = multi_plot(all_results, args.split, plot_mean_step_times_scatter)
        figs.extend(newfigs)


    # Printing
    # ============================================================

    if 'step-times' in args.print_data:

        print_mean_step_times = \
          par(data_print,
              to_print=[('initial_nnode', None),
                        ('-ts', None),
                        ('total_step_time', sp.mean)])

        multi_print(all_results, args.split, print_mean_step_times)


    # End of prints/plots
    # ============================================================


    # if requested then save figures
    # ============================================================

    if args.save_to_dir is not None:

        for fig in figs:

            # Make sure folder exists
            os.makedirs(args.save_to_dir, exist_ok=True)

            name = safefilename(name_fig(fig))

            fig.savefig(pjoin(args.save_to_dir, name + ".pdf"),
                        bbox_inches='tight',
                        pad_inches=0.0,
                        transparent=True)

        if args.save_data_to_dir:
            for data in all_results:
                full_d = data["-outdir"]

                d = pjoin(args.save_to_dir, os.path.split(full_d)[-1])
                os.makedirs(d, exist_ok=True)

                cp(pjoin(full_d, "trace"), pjoin(d, "trace"))
                cp(pjoin(full_d, "info"), pjoin(d, "info"))

    # Show all plots if requested
    if not args.ssh_mode:
        plt.show()


    return 0


def name_fig(fig):
    """Try to come up with a good name for a figure based on labels.
    """

    main_label = [t.get_text() for t in fig.texts if t.get_text() != ""]

    axes_labels = [ax.get_ylabel() + "vs" + ax.get_xlabel()
                   for ax in fig.axes]

    final_label = "-".join(it.chain(main_label, axes_labels))

    return final_label


def safefilename(filename):
    """Strip bad filename characters from string."""
    allowed_symbols = ['_', '-', '.']
    return ''.join(c for c in filename if c.isalnum() or c in allowed_symbols)


# If this script is run from a shell then run main() and return the result.
if __name__ == "__main__":
    sys.exit(main())
