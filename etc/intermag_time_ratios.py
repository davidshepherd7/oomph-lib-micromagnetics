#!/usr/bin/env python3

# Imports from main libraries
import subprocess as subp
import sys
import argparse
import os
import os.path


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


def main():

    # Parse arguments
    # ============================================================

    parser = argparse.ArgumentParser(description=main.__doc__,

    # Don't mess up my formating in the help message
    formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--dir', '-d', action='append',
                        help='Set the directory to look for data in (default "results").')

    parser.add_argument('--split', '-s', action='append',
                        help="Split into different plots for different values of these keys, for keys begining with dash specify them as: `-s='-dt'` to avoid issues with `-` being read as a new argument.")

    parser.add_argument('--skip-failed', action='store_true',
                        help='Skip runs which failed (dir contains file named failed)')


    args = parser.parse_args()

    if (args.dir is None) or (args.dir == []):
        print("No directories given, so parsing ./results")
        args.dir = ["results"]

    if args.split is None:
        args.split = ['mesh', 'h-app', 'initial-m', 'mag-params', 'scale']


    # Parse data
    # ============================================================

    # Get the results that aren't just empty
    really_all_results = mm.parse_parameter_sweep(args.dir,
                                                  skip_failed=args.skip_failed)
    all_results = [d for d in really_all_results if d is not None]


    print(len(all_results), "data sets out of", len(really_all_results), "used",
          "(any others didn't have enough time steps finished).")

    # print("Splitting plots based on values of", args.split)


    # split_data = mm.split_up_stuff(all_results, args.split)

    # Get mean step times
    for data in all_results:
        data['mean_step_time'] = sp.mean(data['total_step_time'])

    # Function specifying which dict keys we care about when splitting into
    # groups for comparison of computation time per step
    def key_filter(k):
        return all([mm.is_arg_key(k),
                    k != "-outdir",
                    k != "-dt",
                    k != "-driver",
                    k != "-fd-jac",
                    k != "-damping",
                    ])


    # Split data up
    grouped = mm.split_to_comparable_groups(all_results, '-ts',
                                            key_filter=key_filter)


    # Print it
    results = []
    for group in grouped:
        # for g in group:
        #     print(g['-ts'],g ['mean_step_time'])
        # print()

        rk_time = sp.mean([d['mean_step_time'] for d in group if d['-ts'] == 'rk2'])
        imr_time = sp.mean([d['mean_step_time'] for d in group if d['-ts'] == 'midpoint-bdf'])

        r = {
             '-ms-method' : group[0]['-ms-method'],
             '-ref' : group[0]['-ref'],
             'ratio' : imr_time/rk_time,
             }

        results.append(r)

    # pprint(results)

    final1 = sp.mean([d['ratio'] for d in results if d['-ms-method'] == 'decoupled'])
    final2 = sp.mean([d['ratio'] for d in results if d['-ms-method'] == 'disabled'])

    # Ignore timing results for 'sphere' because I never got around to
    # writing the code to include the sphere's hms in the Jacobian, so it
    # ends up taking a large number of newton steps to converge.

    # Additional note: some time per step is spent on output. The "full"
    # output (doc_solution()) is only run after certain amounts of
    # simulated time, so the slowdown caused by it is equivalent for
    # implicit/explicit methods. The "partial" output however is called
    # every step, however using valgrind --tool=callgrind we can see that
    # this only accounts for ~5% of execution time, so doesn't really
    # matter for the purposes of our experiment,

    print("decoupled ms time per step ratio:", final1)
    print("disabled ms time per step ratio:", final2)


    return 0


if __name__ == "__main__":
    sys.exit(main())
