#!/bin/bash

# Script to build, run and plot results

# move old results
mv RESLT/* RESLT_old

# The program name - by default this is the same as the folder name with _driver appended
# $1:-.... means if $1 (the first input argument) is set use that, otherwise use ...
# {PWD##/*} gets the present working directory and removes everything before the final / (i.e. just leaving current folder name)
$NAME=${1:-{PWD##/*}_driver}

# make and run other code if make was succesful
if make
then
    
    # Run program and save output to trace file
    ./{$NAME} > RESLT/{$NAME}_trace

    # # Print a solution file to stdout
    # cat RESLT/soln0.dat

    # # Convert to .vtu format (using oomph-convert python script)
    # oomph-convert RESLT/soln* > RESLT/oomph-convert_trace

    # # Open with Paraview
    # paraview --data="RESLT/soln..vtu"

    # zero pad files (to increase length of zero padding change the 3 in "%03d") - this allows gnuplot to plot them in the correct order
    rename 's/\d+/sprintf("%03d",$&)/e' RESLT/soln*

    # plot results using plotting script
    ./plot-gif.sh RESLT/soln*.dat
    
else
    # If make failed then announce and exit
    echo "Build Failed"
    exit 1
fi
