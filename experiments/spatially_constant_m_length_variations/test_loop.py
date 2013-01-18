#!/usr/bin/env python


import subprocess as subp
from multiprocessing import Pool
import itertools
import sys

program='./spatially_constant_m_length_variations'

def dirname_from_parameters(parameters):
    """Convert a list of parameters to a bash-safe(ish) directory name."""

    base = "loop_results/sp"
    for p in parameters:
         base += "_" + str(p)

    # Replace characters that could go wrong in bash/c, probably lots
    # more needed here...
    return base.replace('-','m')

def myrun(parameter_list_as_strings):
    """Make results directory, run the program with the given parameters
    and store the trace.
    """

    outdir = dirname_from_parameters(parameter_list_as_strings)

    subp.call(['mkdir', '-p', outdir])
    tracefile = open(outdir+"/trace", 'w')

    # Write parameter list to file
    paramfile = open(outdir+'/parameters', 'w')
    for p in parameter_list_as_strings:
        paramfile.write(p + " ")
    paramfile.write('\n')
    paramfile.close()

    [eps, dt, mconstraintmethod, damping, hk] = parameter_list_as_strings

    print "Running",eps,dt,mconstraintmethod,damping,hk

    subp.call([program,
             '-eps', eps,
             '-mconstraintmethod', mconstraintmethod,
             '-damp', damping,
             '-hk', hk,
             '-outdir', outdir,
             '-dt', dt,
             '-tmax', '1.0001'],
             stdout=tracefile)

    tracefile.close()

    return


def main():

    # Call make, just in case
    subp.call(["make"])

    # Input lists of parameters to be used
    epsilons = [0]
    dts = [0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001]
    mconstraintmethods = [0, 1, 2]
    dampings = [1.0, 0.5, 0.1, 0.05, 0.01, 0.005]
    hks = [0.0]

    # Convert to a list of lists of strings ready for using call and
    # map_async.
    arguments = [epsilons, dts, mconstraintmethods, dampings, hks]
    args_as_strings = map(lambda x: map(str, x), arguments)

    # Create a list (well, actually an iterator but it doesn't matter) of
    # all possible combinations of arguments (i.e. equivalent to nesting
    # for loops over all argument lists).
    it = itertools.product(*args_as_strings)

    # Start worker processes then run function "myrun" using them
    pool = Pool(processes=7)
    pool.map_async(myrun, it)   # Run all instances
    pool.close()                # Tell python there are no more jobs coming
    pool.join()                 # Wait for everything to finish

    return 0 # Exit successfully

if __name__ == "__main__":
    sys.exit(main())
