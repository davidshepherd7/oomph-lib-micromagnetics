#!/bin/bash

# Script to build, run and plot results

# make and run code if make was succesful
if make
then

    ./demag_driver # > runinfo


# # Combine solutions into single files
#     cat soln[0-9].dat > soln.dat
#     cat exact_soln[0-9].dat > exact_soln.dat
#     cat error[0-9].dat > error.dat
#     rm soln[0-9].dat exact_soln[0-9].dat error[0-9].dat # Delete old files


# Plot results using gnuplot script
    ./plot.sh
    
else
    # If make failed then announce and exit
    echo "Build Failed"
    exit 1
fi
