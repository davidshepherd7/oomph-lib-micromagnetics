#!/bin/bash

# Script to build, run and plot results

# move old results
mv RESLT/* RESLT_old

# make and run other code if make was succesful
if make
then
    
    # Run program
    ./one_d_micromag_driver > one_d_micromag_driver_trace

    # # Print a solution file
    # cat RESLT/soln0.dat

    # # Convert to .vtu format (using oomph-convert python script)
    # oomph-convert RESLT/soln* > RESLT/oomph-convert_trace

    # # Open with Paraview
    # paraview --data="RESLT/soln..vtu"

    # zero pad files (to increase length of zero padding change the 3 in "%03d"
    rename 's/\d+/sprintf("%03d",$&)/e' RESLT/soln*

    # plot results using plotting script
    ./plot-gif.sh RESLT/soln*.dat
    
else
    # If make failed then announce and exit
    echo "Build Failed"
    exit 1
fi
