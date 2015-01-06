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
import glob
import shutil
import multiprocessing as mp
import shlex

# Imports for specific functions
from functools import partial as par
from os.path import join as pjoin
from pprint import pprint
from shutil import copy as cp

# Make sure *this* versions oomphpy is in the path (before any other
# versions in other places)
sys.path.insert(1, pjoin(os.path.dirname(__file__), "../etc"))
import oomphpy
import oomphpy.utils as utils
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


def axis_label_thesisify(label):
    """Replace labels with TeX versions for thesis.
    """
    thesis_labels = {'times': r'$t$',
                     'lnodal' : 'nodal quadrature',
                     'gauss' : 'Gaussian quadrature',
                     'dts' : r'$\Delta_n$',
                     'm length error maxes' : r'$\mathcal{E}_{|\mathbf{m}|}$',
                     '-tol' : r'$\epsilon_\Delta$',
                     'trace values' : r'$y_n$',
                     'error norms' : r'$|y(t_n) - y_n|$',
                     'initial nnode' : r'$N$',

                     'mean of "n solver iters"' : r'mean iterations',
                     'mean of "real solver times"' : r'mean solve time (s)',

                     # newton res stuff
                     'max of "m length error maxes"' : r'$\max (\mathcal{E}_{|\mathbf{m}|})$',
                     'max min of "newton residuals"' : r'max$(||\mathbf{r}_{\mathrm{final}}||_\infty)$',
                     'mean min of "newton residuals"' : r'mean$(||\mathbf{r}_{\mathrm{final}}||_\infty)$',

                     # ode aimr stuff
                     'mean of "dts"' : r'mean($\Delta_n$)',
                     'max of "error norms"' : r'$\max_n |y(t_n) - y_n|$',
                     'maxabs of "energy change"' : r'max $|E_0 - E_n|$',

                     # magnetisation stuff
                     'mean mxs' : r'mean($m_x$)',
                     'mean mys' : r'mean($m_y$)',
                     'mean mzs' : r'mean($m_z$)',

                     'get2 of "trace values"' : r'$m_x(\mathbf{x} = \mathbf{0})$',
                     'get3 of "trace values"' : r'$m_y(\mathbf{x} = \mathbf{0})$',
                     'get4 of "trace values"' : r'$m_z(\mathbf{x} = \mathbf{0})$',

                     # integrals of error norms
                     'aux err1 integral' : r'$\int \, \mathcal{E}_{m_z} \, \mathrm{d}t$',
                     'aux err0 integral' : r'$\int \, \mathcal{E}_{p} \, \mathrm{d}t$',
                     'error norm integral' : r'$\int \, \mathcal{E}_{\mathbf{m}} \, \mathrm{d}t$',
                     'first of "error norms"' : r'$\mathcal{E}_{\mathbf{m}, 1}$',
                     'fakemean of "dts"' : r'$\Delta_n$',
                     }

    # If it matches then change it
    for k, v in thesis_labels.items():
        if k == label:
            return v

    # Otherwise no match: don't change
    return label

def legend_thesisify(label):
    """Replace legends with nice versions for thesis.

    This is an awful hack, if e.g. tr is inside another word the whole
    thing will get replaced... If e.g. tr is inside a replacement the
    result is non-determinisitc, depends on order of dict...  :(
    """

    replaces = {'imr' : 'IMR',
                'tr' : 'TR',
                'bdf2' : 'BDF2',
                'bdf1' : 'BDF1',
                'som-main-block-drop-p ' : r'$\mathcal{P}_2$ ',
                'som-main-block ' : r'$\mathcal{P}_3$ ',
                'gauss' : 'Gaussian',
                'lnodal' : 'nodal',

                # dodgy:
                '$x$' : '$\Delta_n$',
                '$x^2$' : '$\Delta_n^2$',
                '$x^3$' : '$\Delta_n^3$',

                # meshes
                'sq_square' : 'square',
                'st_square' : 'triangular',
                }

    for k, v in replaces.items():
        if k in label:
            label = label.replace(k, v)
    else:
        return label



def latex_safe(string):
    """Strip out characters that will interfere with latex.
    * _ -> space
    * ??
    """
    return string.replace("_", " ")


def make_label(p, op=None, op_name=None):
    """Create a label from field name and operation.
    """

    # Use the most specific name we have
    if op_name is not None:
        return op_name + " of " + '"' + p + '"'
    elif op is not None:
        return op.__name__ + " of " + '"' + p + '"'
    else:
        return p


def make_axis_label(*args, **kwargs):
    return latex_safe(make_label(*args, **kwargs))


def make_symbol_iter(symbol_type):
    if symbol_type == "symbols":
        symbols = ['x', 'o', '+', '^', '*', 's', ]
    elif symbol_type == "lines":
        symbols = ['-', '--', '-.', ':']
    else:
        raise ValueError

    return it.cycle(symbols)

def make_colour_iter():
    colours = ['r', 'g', 'b', 'k', 'c', 'm', 'y']
    return it.cycle(colours)


# Plotting functions
# ============================================================


def plot_vs_thing(xthing, data, plot_values,
                  operations_on_values=None,
                  operation_on_x=None,
                  labels=None,
                  y_axis_lims=None,
                  xscale="linear",
                  yscale="linear",
                  skip_first_n=0,
                  plot_points=False):
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
    if all([d['-tol'] == 0 or d['-dummy-adaptivity'] == 1 for d in data]):
        dt_label = '-dt'
    else:
        dt_label = '-tol'

    # Default labels
    if labels is None:
        labels = ['-ref', '-ts', dt_label]

    for axes, p, op in zip(axesarray, plot_values, operations_on_values):

        if plot_points:
            symbols = make_symbol_iter("symbols")
            lines = it.repeat('-')
        else:
            lines = make_symbol_iter("lines")
            symbols = it.repeat(None)
        colours = make_colour_iter()

        for d in data:

            name = " ".join([str(d[l]) for l in labels])

            initial_vals = d.get(p)
            if initial_vals is not None:

                # Do things to y-values (e.g. mean, max, ...)
                if op is not None:
                    vals = list(map(op, d[p]))

                else:
                    vals = d[p]

                # Do things to x-values (e.g. convert units)
                if operation_on_x is not None:
                    xs = list(map(operation_on_x, d[xthing]))
                else:
                    xs = d[xthing]

                # Skip first few steps if requested
                xs = xs[skip_first_n:-1]
                vals = vals[skip_first_n:-1]

                colour = next(colours)

                axes.plot(xs, vals,
                          linestyle=next(lines),
                          marker=next(symbols),
                          color=colour,
                          markerfacecolor='none', #
                          markeredgecolor=colour,
                          label=name)
            else:
                sys.stderr.write("Not plotting " + p
                                 + " because I couldn't find the data needed.\n")

                # shift to next symbol + colour anyway for consistency
                next(symbols)
                next(lines)
                next(colours)

        axes.set_ylabel(make_axis_label(p, op))

        axes.set_xscale(xscale)
        axes.set_yscale(yscale)

    # x label only on last axis
    axesarray[-1].set_xlabel(xthing)

    # Add a single legend
    axesarray[0].legend(loc=0)

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
               xscale="log",
               yscale="log",
               y_operation_name=None,
               x_operation_name=None,
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


    symbols = make_symbol_iter("symbols")
    colours = make_colour_iter()

    split_data = utils.split_up_stuff(data, dataset_split_keys)


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

        colour = next(colours)

        # Plot
        ax.plot(xs, ys,
                marker=next(symbols), alpha=0,
                markerfacecolor='none',
                markeredgecolor=colour,
                markersize=8,
                label=name)

    # Label
    ax.set_xlabel(make_axis_label(x_value, x_operation, x_operation_name))
    ax.set_ylabel(make_axis_label(y_value, y_operation, y_operation_name))

    ax.set_xscale(xscale)
    ax.set_yscale(yscale)


    # Pick the right axis scales
    ax.axis('auto')

    ax.legend(loc=0)



    return fig


def multi_plot(data, keys_to_split_on, plot_function):
    """Split data into sub-sets (based on keys_to_split_on) and plot each
    subset as specified by (partially evaluated) function given.
    """

    # Divide into separate data sets
    split_data = utils.split_up_stuff(data, keys_to_split_on)

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
    split_data = utils.split_up_stuff(data, keys_to_split_on)

    for dataset in split_data:

        # Make a title based on the keys which specify this data set (all
        # data in dataset have the same value for the keys in
        # keys_to_split_on so we just get it from the first one).
        labels = []
        for k in keys_to_split_on:
            try:
                this_str = str(dataset[0][k])
                labels.append(str(k) + " = " + this_str + "; ")

            # Ignore keys that don't exist
            except KeyError:
                pass

        print(''.join(labels))

        # Print data
        print_function(dataset)


    return


def data_print(datasets, to_print, delim="; ", labels=None):

    if labels is None:
        labels = []

    # Also print labels
    to_print = to_print + [(l, None) for l in labels]

    # Print headers for each column
    for a, op in to_print:
        print(make_label(a, op), end=delim)

    print()

    # Print the values for each dataset
    for d in datasets:
        for a, op in to_print:
            if op is None:
                print(d[a], end=delim)
            else:
                print(op(d[a]), end=delim)

        print()

    print()


def shift_relaxation_times(data):
    """Modify times in data such that they are continuous even
    when using relax(..).
    """

    ts = data['times']
    dts = data['dts']

    second_t0_i = None
    for i, (t, dt) in enumerate(zip(ts[3:], dts[3:]), 3):

        # Look for zeros
        if abs(t) < 1e-10:
            second_t0_i = i
            t_relax_max = ts[i-1]
            break

    # if we found some t = 0 not in the 0th place in the list
    if second_t0_i is not None:

        # Make the relaxation t values negative
        new_ts = sp.array([t - t_relax_max for t in ts[0:second_t0_i]] + list(ts[second_t0_i:]))

        # Check it did what I think it did...
        assert abs(new_ts[second_t0_i]) < 1e-10

        data['times'] = new_ts

    return


def maxabs(l):
    return max((abs(x) for x in l))

def datasort(data1, data2):
    """Define an ordering on "data" dictionaries.

    Based on the -ts argument then the -renormalise argument, this
    should allow consistent colours in plots.
    """

    return


def main():
    """

    """

    # Parse inputs
    # ============================================================

    parser = argparse.ArgumentParser(description=main.__doc__,
                                     # Show defaults in help
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('dir', nargs="*", default=["results"],
                        help='Set the directory(s) to look for data in')

    parser.add_argument('--print-all-data', action='store_true',
                        help='Pretty print all data to stdout')

    parser.add_argument('--plots', '-p', action='append', default=[],
                        help='choose what to plot')

    parser.add_argument('--print-data', '-v', action='append', default=[],
                        help='Choose values to print')

    parser.add_argument('--label', '-l', action='append',
                        help='Add additional labels to line, for keys begining with dash specify them as: `-s=-dt` to avoid issues with `-` being read as a new argument.')

    parser.add_argument('--split', '-s', action='append', default=[],
                        help="Split into different plots for different values of these keys, for keys begining with dash specify them as: `-s=-dt` to avoid issues with `-` being read as a new argument.")

    parser.add_argument('--scatter-split', '-t', action='append',
                        help="Split into different data sets in a scatter plot for different values of these keys, for keys begining with dash specify them as: `-s=-dt` to avoid issues with `-` being read as a new argument.")

    parser.add_argument('--skip-failed', action='store_true',
                        help='Skip runs which failed (dir contains file named failed)')

    parser.add_argument('--filter', '-f', action='append', default=[],
                        help="""Filter runs based on parameters. Input should be a python pair like ('-tol', 1e-3) where the first entry is the dictionary key to check and the second entry is the value it should take. Note the quotes around the key.""")

    parser.add_argument('--not-filter', action='append', default=[],
                        help="""Oposite of -filter: only use runs with other values of this parameter. Input should be a python pair like ('-tol', 1e-3) where the first entry is the dictionary key to check and the second entry is the value it should not take. Note the quotes around the key.""")

    parser.add_argument('--save-to-dir',
                        help='Save figures as pdfs into the specified folder.')

    parser.add_argument('--ssh-mode', action='store_true',
                        help="Vastly speed up operation over ssh: don't show plots and use Agg backend.")

    parser.add_argument('--save-data-to-dir', action='store_true',
                        help='Also copy trace and info files to the directory.')

    parser.add_argument('--thesis', action='store_true',
                        help='Do some modifications to make plots prettier.')

    parser.add_argument('--debug-mode', action='store_true',
                        help = 'Enable debugging mode (run in serial).')

    parser.add_argument('--ncores', '-n', '-j', default=mp.cpu_count(), type=int,
                        help='Number of cores to use')

    parser.add_argument('--points', action='store_true',
                        help='Plot points on line graphs')

    args = parser.parse_args()


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
                                                  skip_failed=args.skip_failed,
                                                  serial_mode=args.debug_mode,
                                                  processes=args.ncores)
    all_results = [d for d in really_all_results if d is not None]

    print(len(all_results), "data sets out of", len(really_all_results), "used",
          "(any others didn't have enough time steps finished).")

    print("Splitting plots based on values of", args.split)
    print("Labeling  plots based on values of", args.label)

    # Print if needed
    if args.print_all_data:
        pprint(all_results)

    def literal_eval_hack(input):
        """A literal eval that can handle being given numbers instead of

        a string containing the text of a number.
        """
        try:
            return ast.literal_eval(input)
        except ValueError:
            pass

        try:
            return ast.literal_eval(str(input))
        except ValueError:
            pass

        try:
            return ast.literal_eval('"' + str(input) + '"')
        except ValueError:
            print("failed to eval:", input)
            raise


    # Filter the results based on the arguments given, using literal eval
    # means that floats are compared as floats, strings as strings etc.
    try:
        for f in args.filter:
            key, value = ast.literal_eval(f)
            all_results = [d for d in all_results if literal_eval_hack(str(d[key])) == value]
            print("filtering with", f, ".", len(all_results), "results left")


        for f in args.not_filter:
            key, value = ast.literal_eval(f)
            all_results = [d for d in all_results if literal_eval_hack(str(d[key])) != value]
            print("filtering with not", f, ".", len(all_results), "results left")

    except ValueError:
        print("failed on:", f, key, value, all_results[0][key])
        raise


    for data in all_results:
        shift_relaxation_times(data)


    # Define a plot function including the cli specified options
    def plot_vs_time_with_cli_args(*pargs, **kwargs):
        return plot_vs_time(*pargs,
                            labels=args.label,
                            plot_points=args.points,
                            **kwargs)



    # Do actual plots
    # ============================================================

    # Storage for figures
    figs = []


    # Plot error norm vs time
    if 'err' in args.plots:
        plot_errors = par(plot_vs_time_with_cli_args,
                          plot_values=['error_norms','dts', 'trace_values'])

        newfigs = multi_plot(all_results, args.split, plot_errors)
        figs.extend(newfigs)


    if 'aux_err0' in args.plots:
        plot_errors = par(plot_vs_time_with_cli_args,
                          plot_values=['aux_error_norms','dts', 'trace_values'],
                          operations_on_values=[lambda x:x[0], None, None])

        newfigs = multi_plot(all_results, args.split, plot_errors)
        figs.extend(newfigs)

    if 'aux_err1' in args.plots:
        plot_errors = par(plot_vs_time_with_cli_args,
                          plot_values=['aux_error_norms','dts', 'trace_values'],
                          operations_on_values=[lambda x:x[1], None, None])

        newfigs = multi_plot(all_results, args.split, plot_errors)
        figs.extend(newfigs)

    if 'swerr' in args.plots:
        plot_errors = par(plot_vs_time_with_cli_args,
                          plot_values=['switching_time_error','dts', 'trace_values'])

        newfigs = multi_plot(all_results, args.split, plot_errors)
        figs.extend(newfigs)

    # Plot m averages vs time
    if 'm' in args.plots:
        plot_m_averages = par(plot_vs_time_with_cli_args,
                              plot_values=['mean_mxs','mean_mys','mean_mzs','dts'],
                              y_axis_lims=[[-1,1], [-1,1], [-1,1], None])

        newfigs = multi_plot(all_results, args.split, plot_m_averages)
        figs.extend(newfigs)

    if 'm-field' in args.plots:
        plot_m_averages = par(plot_vs_time_with_cli_args,
                              plot_values=['mean_mxs','mean_mys','mean_mzs','dts',
                                           'h_applied_first_element'],
                              y_axis_lims=[[-1,1], [-1,1], [-1,1], None, None])

        newfigs = multi_plot(all_results, args.split, plot_m_averages)
        figs.extend(newfigs)


    # kind of assumes 3 trace values
    if 'tracem' in args.plots:

        # Named functions so we can replace with good values for thesis version
        def get2(y):
            return y[2]
        def get3(y):
            return y[3]
        def get4(y):
            return y[4]

        plot_traces = par(plot_vs_time_with_cli_args,
                          plot_values=['trace_values','trace_values','trace_values', 'dts'],
                          operations_on_values=[get2, get3, get4, None])

        newfigs = multi_plot(all_results, args.split, plot_traces)

        for fig in newfigs:
            for ax in fig.axes:
                ax.set_xlim(-0.1, None)

        figs.extend(newfigs)


    if 'trace' in args.plots:
        plot_traces = par(plot_vs_time_with_cli_args,
                          plot_values=['trace_values', 'dts'])

        newfigs = multi_plot(all_results, args.split, plot_traces)
        figs.extend(newfigs)

    if 'energy' in args.plots:
        plot_traces = par(plot_vs_time_with_cli_args,
                          plot_values=['total_energy',
                                        'exchange_energy',
                                        'zeeman_energy',
                                        'crystalline_anisotropy_energy',
                                        'magnetostatic_energy'])

        newfigs = multi_plot(all_results, args.split, plot_traces)
        figs.extend(newfigs)

    if 'newt' in args.plots:
        plot_newton_iters = par(plot_vs_time_with_cli_args,
                                plot_values=['n_newton_iters','dts'])

        newfigs = multi_plot(all_results, args.split, plot_newton_iters)
        figs.extend(newfigs)

    if 'allits' in args.plots:
        plot_all_iters = par(plot_vs_time_with_cli_args,
                             plot_values=['n_newton_iters', 'n_solver_iters', 'preconditioner_setup_times', 'dts'],
                             operations_on_values=[None, sp.mean, sp.mean, None])

        newfigs = multi_plot(all_results, args.split, plot_all_iters)
        figs.extend(newfigs)

    if 'soltimes' in args.plots:
        plot_sol_time_averages = par(plot_vs_time_with_cli_args,
                                plot_values=['solver_times','jacobian_setup_times'],
                                operations_on_values=[sp.mean, sp.mean])

        newfigs = multi_plot(all_results, args.split, plot_sol_time_averages)
        figs.extend(newfigs)

    # Plot |m| error vs time
    if 'ml' in args.plots:
        plot_ml_error_vs_time = par(plot_vs_time_with_cli_args,
                                    plot_values=['m_length_error_maxes', 'dts'])

        newfigs = multi_plot(all_results, args.split, plot_ml_error_vs_time)
        figs.extend(newfigs)

        # Plot |m| error vs time
    if 'ml-only' in args.plots:
        plot_ml_error_vs_time = par(plot_vs_time_with_cli_args,
                                    plot_values=['m_length_error_maxes'])

        newfigs = multi_plot(all_results, args.split, plot_ml_error_vs_time)
        figs.extend(newfigs)

    if 'ml-only-log' in args.plots:
        plot_ml_error_vs_time = par(plot_vs_time_with_cli_args,
                                    plot_values=['m_length_error_maxes'],
                                    yscale="log")

        newfigs = multi_plot(all_results, args.split, plot_ml_error_vs_time)
        figs.extend(newfigs)

    if 'energy-only-log' in args.plots:
        plot_ml_error_vs_time = par(plot_vs_time_with_cli_args,
                                    plot_values=['energy_change'],
                                    operations_on_values=[abs],
                                    yscale="log")

        newfigs = multi_plot(all_results, args.split, plot_ml_error_vs_time)
        figs.extend(newfigs)

    if 'mumag4-energy-only-log' in args.plots:
        plot_ml_error_vs_time = par(plot_vs_time_with_cli_args,
                                    plot_values=['energy_change'],
                                    operations_on_values=[abs],
                                    yscale="log")

        newfigs = multi_plot(all_results, args.split, plot_ml_error_vs_time)

        for fig in newfigs:
            for ax in fig.axes:
                ax.set_xlim(0, None)

        figs.extend(newfigs)

    if 'mumag4-comparison-plot' in args.plots:

        def to_nanoseconds(t):
            return t * (5.653e-3)

        plot_ml_error_vs_time = par(plot_vs_time_with_cli_args,
                                    plot_values=['mean_mys'],
                                    operation_on_x=to_nanoseconds)

        newfigs = multi_plot(all_results, args.split, plot_ml_error_vs_time)

        for fig in newfigs:
            for ax in fig.axes:
                ax.set_xlim(0, 1.0)
                ax.set_ylim(-0.5, 0.5)


        figs.extend(newfigs)

    if 'lte' in args.plots:
        plot_ml_error_vs_time = par(plot_vs_time_with_cli_args,
                                    plot_values=['LTE_norms', 'dts'],
                                    skip_first_n=1)
        newfigs = multi_plot(all_results, args.split, plot_ml_error_vs_time)
        figs.extend(newfigs)


    # Plot solver iterations vs steps
    if 'its' in args.plots:
        plot_iters_step = par(plot_vs_time_with_cli_args, plot_values=['dts', 'n_solver_iters'],
                              operations_on_values=[identity, sp.mean])

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
                                  y_operation=damping_error_mean)
        newfigs = multi_plot(all_results, args.split, plot_damping_errors)
        figs.extend(newfigs)


    # Plot error in effective damping vs step size
    if 'damping' in args.plots:

        plot_damping_errors = par(plot_vs_time_with_cli_args,
                                  plot_values=[# 'rel_damping_error',
                                               # 'abs_damping_error',
                                               'alt_effective_damping',
                                               'dts',
                                               ],
                                  skip_first_n=1, # first data points are wrong
                                  )

        newfigs = multi_plot(all_results, args.split, plot_damping_errors)
        figs.extend(newfigs)


    if 'wc-time' in args.plots:
        plot_wall_time_vs_time = \
          par(plot_vs_time_with_cli_args,
              plot_values=['unix_timestamp_diffs', 'dts'])

        newfigs = multi_plot(all_results, args.split, plot_wall_time_vs_time)
        figs.extend(newfigs)

    if 'scatter-step-time' in args.plots:
        step_time_scatter = \
          par(my_scatter,
              labels=args.label,
              dataset_split_keys=args.scatter_split,
              y_value='total_step_time',
              yscale='linear',
              x_value='initial_nnode',
              y_operation=sp.mean)


        newfigs = multi_plot(all_results, args.split, step_time_scatter)
        figs.extend(newfigs)


    if 'scatter-aux_err0-integrals' in args.plots:
        plot_errors = par(my_scatter,
                          labels=args.label,
                          dataset_split_keys=args.scatter_split,
                          x_value='dts',
                          y_value='aux_err0_integral',
                          x_operation=sp.mean,
                          x_operation_name="fakemean")

        newfigs = multi_plot(all_results, args.split, plot_errors)


        figs.extend(newfigs)

    if 'scatter-aux_err1-integrals' in args.plots:
        plot_errors = par(my_scatter,
                          labels=args.label,
                          dataset_split_keys=args.scatter_split,
                          x_value='dts',
                          y_value='aux_err1_integral',
                          x_operation=sp.mean,
                          x_operation_name="fakemean")

        newfigs = multi_plot(all_results, args.split, plot_errors)

        for fig in newfigs:
            for ax in fig.axes:
                xs = sp.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 5)
                ax.plot(xs, [x**2 for x in xs], 'k-', label="$x^2$")
                ax.legend(loc=0)

        figs.extend(newfigs)

    if 'scatter-error-integral' in args.plots:
        plot_errors = par(my_scatter,
                          labels=args.label,
                          dataset_split_keys=args.scatter_split,
                          x_value='dts',
                          y_value='error_norm_integral',
                          x_operation=sp.mean,
                          x_operation_name="fakemean")

        newfigs = multi_plot(all_results, args.split, plot_errors)

        figs.extend(newfigs)

    if 'scatter-error-period' in args.plots:

        def thingy(values):
            phase_err = sp.array([m[0] - 2*sp.pi for m in values])
            inds, zero_crossings = mm.find_zero_crossings(phase_err)

            return len(phase_err)*4*sp.pi + phase_err[-1] + 2*sp.pi

        plot_errors = par(my_scatter,
                          labels=args.label,
                          dataset_split_keys=args.scatter_split,
                          x_value='dts',
                          y_value='aux_error_norms',
                          y_operation=thingy,
                          x_operation=sp.mean,
                          x_operation_name="fakemean")

        newfigs = multi_plot(all_results, args.split, plot_errors)

        figs.extend(newfigs)


    if 'scatter-step-time-log' in args.plots:
        step_time_scatter = \
          par(my_scatter,
              labels=args.label,
              dataset_split_keys=args.scatter_split,
              y_value='total_step_time',
              yscale='log',
              x_value='initial_nnode',
              y_operation=sp.mean)


        newfigs = multi_plot(all_results, args.split, step_time_scatter)
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

    if 'scatter-mean-solve-times-nnode' in args.plots:
        # for d in all_results:
            # print(sp.mean(d['total_step_time']))

        plot_mean_step_times_scatter = \
          par(my_scatter,
              labels=args.label,
              dataset_split_keys=args.scatter_split,
              yscale='linear',
              x_value='initial_nnode',
              y_value='real_solver_times',
              y_operation=lambda x:sp.mean(list(it.chain(*x))),
              y_operation_name="mean")

        newfigs = multi_plot(all_results, args.split, plot_mean_step_times_scatter)
        figs.extend(newfigs)

    if 'scatter-mean-solve-times-dt' in args.plots:
        # for d in all_results:
            # print(sp.mean(d['total_step_time']))

        plot_mean_step_times_scatter = \
          par(my_scatter,
              labels=args.label,
              dataset_split_keys=args.scatter_split,
              yscale='linear',
              x_value='dts',
              y_value='real_solver_times',
              x_operation=sp.mean,
              y_operation=lambda x:sp.mean(list(it.chain(*x))),
              y_operation_name="mean")

        newfigs = multi_plot(all_results, args.split, plot_mean_step_times_scatter)
        figs.extend(newfigs)

    if 'scatter-mean-prec-times-nnode' in args.plots:
        # for d in all_results:
            # print(sp.mean(d['total_step_time']))

        plot_mean_step_times_scatter = \
          par(my_scatter,
              labels=args.label,
              dataset_split_keys=args.scatter_split,
              yscale='linear',
              x_value='initial_nnode',
              y_value='preconditioner_setup_times',
              y_operation=lambda x:sp.mean(list(it.chain(*x))),
              y_operation_name="mean")

        newfigs = multi_plot(all_results, args.split, plot_mean_step_times_scatter)
        figs.extend(newfigs)

    if 'scatter-its-dt' in args.plots:
        plot_mean_step_times_scatter = \
          par(my_scatter,
              labels=args.label,
              dataset_split_keys=args.scatter_split,
              yscale='linear',
              x_value='dts',
              y_value='n_solver_iters',
              x_operation=sp.mean,
              y_operation=lambda x:sp.mean(list(it.chain(*x))),
              y_operation_name="mean")

        newfigs = multi_plot(all_results, args.split, plot_mean_step_times_scatter)
        figs.extend(newfigs)

        # iters axis should start from 0
        for fig in newfigs:
            for ax in fig.axes:
                ax.set_ylim(bottom=0)

    if 'scatter-its' in args.plots:
        plot_mean_step_times_scatter = \
          par(my_scatter,
              labels=args.label,
              dataset_split_keys=args.scatter_split,
              yscale='linear',
              x_value='initial_nnode',
              y_value='n_solver_iters',
              y_operation=lambda x:sp.mean(list(it.chain(*x))),
              y_operation_name="mean")

        newfigs = multi_plot(all_results, args.split, plot_mean_step_times_scatter)
        figs.extend(newfigs)

        # iters axis should start from 0
        for fig in newfigs:
            for ax in fig.axes:
                ax.set_ylim(bottom=0)

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

        for fig in newfigs:
            for ax in fig.axes:
                xs = sp.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 5)
                ax.plot(xs, [x**2 for x in xs], 'k-', label="$x^2$")
                ax.legend(loc=0)

        figs.extend(newfigs)


    if 'scatter-err-dts-noguides' in args.plots:
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


    if 'scatter-denergy-dt' in args.plots:
        plot_err_scatter = \
          par(my_scatter,
              labels=args.label,
              dataset_split_keys=args.scatter_split,
              x_value='dts',
              y_value='energy_change',
              y_operation=maxabs,
              x_operation=sp.mean)

        newfigs = multi_plot(all_results, args.split, plot_err_scatter)
        figs.extend(newfigs)


    if 'scatter-swerr-dts' in args.plots:
        plot_err_scatter = \
          par(my_scatter,
              labels=args.label,
              dataset_split_keys=args.scatter_split,
              x_value='dts',
              y_value='switching_time_error',
              y_operation=max,
              x_operation=sp.mean)

        newfigs = multi_plot(all_results, args.split, plot_err_scatter)
        figs.extend(newfigs)

        for fig in newfigs:
            for ax in fig.axes:
                xs = sp.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 5)
                ax.plot(xs, [x**2 for x in xs], 'k-', label="$x^2$")
                ax.legend(loc=0)


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
              x_value='dts',
              y_value='error_norms',
              y_operation=lambda y:y[1],
              y_operation_name="first",
              x_operation=sp.mean,
              x_operation_name="fakemean")

        newfigs = multi_plot(all_results, args.split, plot_err_scatter)

        for fig in newfigs:
            for ax in fig.axes:
                xs = sp.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 5)
                ax.plot(xs, [100*x**2 for x in xs], 'k-', label="$100x^2$")
                ax.legend(loc=0)

        figs.extend(newfigs)

    if 'scatter-ml-newtres' in args.plots:

        def mean_min(dataset):
            mins = []
            for residuals in dataset:
                if len(residuals) > 0:
                    mins.append(min(residuals))

            return sp.mean(mins)

        fplot = \
          par(my_scatter,
              labels=args.label,
              dataset_split_keys=args.scatter_split,
              y_value='m_length_error_maxes',
              y_operation=max,
              x_value='newton_residuals',
              x_operation=mean_min)

        newfigs = multi_plot(all_results, args.split, fplot)
        figs.extend(newfigs)

    if 'scatter-energy-change-newtres' in args.plots:

        def mean_min(dataset):
            mins = []
            for residuals in dataset:
                if len(residuals) > 0:
                    mins.append(min(residuals))

            return sp.mean(mins)

        fplot = \
          par(my_scatter,
              labels=args.label,
              dataset_split_keys=args.scatter_split,
              y_value='energy_change',
              y_operation=max,
              x_value='newton_residuals',
              x_operation=mean_min)

        newfigs = multi_plot(all_results, args.split, fplot)
        figs.extend(newfigs)


    if 'scatter-ml-maxminnewtres' in args.plots:

        def max_min(dataset):
            mins = []
            for residuals in dataset:
                if len(residuals) > 0:
                    mins.append(min(residuals))

            return max(mins)

        fplot = \
          par(my_scatter,
              labels=args.label,
              dataset_split_keys=args.scatter_split,
              y_value='m_length_error_maxes',
              y_operation=max,
              x_value='newton_residuals',
              x_operation=max_min)

        newfigs = multi_plot(all_results, args.split, fplot)
        figs.extend(newfigs)

    if 'scatter-ml-newttol' in args.plots:

        fplot = \
          par(my_scatter,
              labels=args.label,
              dataset_split_keys=args.scatter_split,
              y_value='m_length_error_maxes',
              y_operation=max,
              x_value='-newton-tol')

        newfigs = multi_plot(all_results, args.split, fplot)
        figs.extend(newfigs)

    if 'scatter-ml-initialnewtres' in args.plots:

        def mean_initial(dataset):
            initial = []
            for residuals in dataset:
                if len(residuals) > 0:
                    initial.append(residuals[0])

            return sp.mean(initial)

        fplot = \
          par(my_scatter,
              labels=args.label,
              dataset_split_keys=args.scatter_split,
              y_value='m_length_error_maxes',
              y_operation=max,
              x_value='newton_residuals',
              x_operation=mean_initial)

        newfigs = multi_plot(all_results, args.split, fplot)
        figs.extend(newfigs)

    if 'scatter-newt' in args.plots:
        plot_mean_step_times_scatter = \
          par(my_scatter,
              labels=args.label,
              dataset_split_keys=args.scatter_split,
              x_value='initial_nnode',
              y_value='n_newton_iters',
              yscale='linear',
              y_operation=sp.mean)

        newfigs = multi_plot(all_results, args.split, plot_mean_step_times_scatter)
        figs.extend(newfigs)


    if 'convergence-c' in args.plots:

        # First we need to add convergence rate to the datasets.
        conv_datasets = mm.split_to_comparable_groups(all_results, '-ref')
        for dataset in conv_datasets:
            rate = mm.convergence_rate(dataset, error_norm_norm=max)

            # Add to all data in the set
            for d in dataset:
                d['convergence_rate'] = rate


        fplot = par(my_scatter,
                    labels=args.label,
                    dataset_split_keys=args.scatter_split,
                    x_value='-wave-solution-c',
                    y_value='convergence_rate',
                    xscale='linear', yscale='linear')

        # Only plot the datasets with the largest refinement (others are
        # equivalent except for refinement level).
        max_ref = max([d['-ref'] for d in all_results])
        max_ref_results = [d for d in all_results if d['-ref'] == max_ref]
        newfigs = multi_plot(max_ref_results, args.split, fplot)



    # Printing
    # ============================================================

    if 'step-times' in args.print_data:

        print_mean_step_times = \
          par(data_print,
              to_print=[('initial_nnode', None),
                        ('-ts', None),
                        ('total_step_time', sp.mean)],
              labels=args.label)


        multi_print(all_results, args.split, print_mean_step_times)


    if 'ml' in args.print_data:

        print_mean_step_times = \
          par(data_print,
              to_print=[('-ts', None),
                        ('m_length_error_maxes', max)],
              labels=args.label)


        multi_print(all_results, args.split, print_mean_step_times)

    # print max(max(m_length_error_maxes)) (over all time, over all
    # parameter sets)
    if 'max-max-ml' in args.print_data:

        def print_max_max(dataset):
            max_max = max([max(d['m_length_error_maxes']) for d in dataset])
            print(max_max)
            return

        multi_print(all_results, args.split, print_max_max)



    # End of prints/plots
    # ============================================================

    # Generate filenames for each figure, and add as a new class
    # member. Must be done before "thesisify" because that changes the axis
    # labels to latex code.
    for fig in figs:
        fig._ds_name = name_fig(fig)


    # Prettyfy for thesis
    if args.thesis:
        for f in figs:

            # No titles on plots
            f.suptitle('')

            # Replace strings in axis labels
            for ax in f.axes:
                ax.set_xlabel(axis_label_thesisify(ax.get_xlabel()))
                ax.set_ylabel(axis_label_thesisify(ax.get_ylabel()))

            # If there is only one line then remove the legend
            for ax in f.axes:
                leg = ax.get_legend()
                if (leg is not None) and (len(leg.get_lines()) == 1):
                    ax.legend_ = None

            # Replace strings in legends
            for ax in f.axes:
                if ax.get_legend() is not None: # if legend exists
                    handles, labels = ax.get_legend_handles_labels()
                    labels = [legend_thesisify(l) for l in labels]
                    ax.legend(handles, labels, loc=0)



    # if requested then save figures
    # ============================================================

    if args.save_to_dir is not None:

        # Make sure folder exists
        os.makedirs(args.save_to_dir, exist_ok=True)

        # Write plot command used to output dir
        plot_command_filename = pjoin(args.save_to_dir, "plot_command")
        with open(plot_command_filename, 'w') as plot_command_file:
            plot_command_file.write(' '.join([shlex.quote(a) for a in sys.argv]))
            plot_command_file.write('\n')


        # Copy any parameter sets in data dirs to output dir:

        # First delete existing parameter sets in folder
        for f in glob.glob(pjoin(args.save_to_dir, "parameter_set*")):
            os.remove(f)

        # Find all "parameter_file" files (non-recursively, recusive would
        # have to search *all* data files, slow..)
        param_files = list(filter(os.path.exists, [pjoin(d, "parameter_file") for d in args.dir]))

        # Copy to output dir (with number appended to prevent overwrites)
        for i, d in enumerate(param_files):
            shutil.copyfile(d, pjoin(args.save_to_dir, "parameter_file") + "_" + str(i))


        # Copy the figures
        for fig in figs:

            # construct path
            name = safefilename(fig._ds_name)
            path = pjoin(args.save_to_dir, name + ".pdf")

            # Do it
            fig.savefig(path,
                        bbox_inches='tight',
                        pad_inches=0.0,
                        transparent=True)

            print("saved plot as", path)


        # Copy the data if requested
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


# If this script is run from a shell then run main() and return the result.
if __name__ == "__main__":
    sys.exit(main())
