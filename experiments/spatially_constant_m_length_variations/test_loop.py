#!/usr/bin/env python


import subprocess as sp
import sys

program='./spatially_constant_m_length_variations'

def mydirname(parameters):
    """Convert a list of strings to a bash-safe dir name.
    """
    mystring = "loop_results/sp"

    for p in parameters:
        mystring += "_" + str(p)

    # lots of others...
    # nasty!
    return mystring.replace('-','m')

if __name__ == "__main__":


    sp.call(["make"])
    sp.call(["rm", "loop_results/*", "-rf"])

    eps_list = map(str, [1e-3, 1e-4, 1e-6])
    renormalise_list = map(str, [0,1])
    damping_list = map(str, [1.0, 0.5, 0.1, 0.05])
    hk_list = map(str, [0.0, 0.5, 1.0])

    for eps in eps_list:
        for renormalise_bool in renormalise_list:
            for damping in damping_list:
                for hk in hk_list:

                    outdir = mydirname([eps,renormalise_bool,damping,hk])
                    sp.call(['mkdir', '-p', outdir])

                    sp.call([program,
                             '-eps', eps,
                             '-renormalise', renormalise_bool,
                             '-damp', damping,
                             '-hk', hk,
                             '-outdir', outdir,
                             ])


    sys.call(['./process_data.sh',''])


    sys.exit(0) # Exit successfully
