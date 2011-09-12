#!/bin/bash

# Script to build, run and plot results

# move old results
mv RESLT/* RESLT_old

# make and run other code if make was succesful
if make
then
    
    # Run program
    ./demag_driver > runtrace

    # Convert to .vtu format (using oomph-convert python script)
    oomph-convert RESLT/soln*

    # Open with Paraview
    paraview --data="RESLT/soln..vtu"
    
else
    # If make failed then announce and exit
    echo "Build Failed"
    exit 1
fi
