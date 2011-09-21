#!/bin/bash

# Script to build, run and plot results

# If any command fails then exit the script
set -o errexit

# move old results
mv RESLT/* RESLT_old

# The program name - by default this is the same as the folder name with _driver appended so get the name of the folder we are in:
NAME=`basename $PWD`_driver

# Build the code
make

# Run program and save output to trace file
./${NAME} > RESLT/${NAME}_trace

# # Print a solution file to stdout
# print 0

# # Convert to .vtu format (using oomph-convert python script)
# convert-vtu

# # Open with Paraview
# paraview --data="RESLT/soln..vtu"

# zero pad files (to increase length of zero padding change the 3 in "%03d") - this allows gnuplot to plot them in the correct order
rename 's/\d+/sprintf("%03d",$&)/e' RESLT/soln*

# plot results
plot-animation RESLT/soln*.dat

################################################################################
# Function definitions
################################################################################

# Print the solution to std out (i.e. the terminal normally)
function print
{
    cat RESLT/soln${1}.dat
}

# Convert all soln*.dat files in RESLT folder into .vtu files ready for use with paraview
# Put all output from the conversion into file oomph-convert_trace
function convert-vtu
{
    oomph-convert RESLT/soln* > RESLT/oomph-convert_trace
}


# Plot an animated gif of the time evolution from the data files, one frame per data file
# ??ds need to set the titles and columns to be used manually
########################################################################
function plot-animation
{
    # Set options
    echo "clear" > .plot.gpi
    echo "reset" >> .plot.gpi
    echo "set terminal gif animate delay 8" >> .plot.gpi # delay = time between frames
    echo "set output \"RESLT/animate.gif\"" >> .plot.gpi # ouput file name

    # set ranges of plots (to ensure they are the same for all frames)
    #echo "set xrange [0:1]" >> .plot.gpi
    echo "set yrange [-1:1]" >> .plot.gpi

    # plot data from each file
    for file in $@
    do
	echo "plot \"$file\" using 1:2 t \"\phi\", \"$file\" using 1:3 title \"M_x\", \"$file\" using 1:4 t \"M_y\", \"$file\" using 1:5 t \"M_z\" " >> .plot.gpi
    done

    # Run the commands
    gnuplot .plot.gpi

    # Show the gif
    eog RESLT/animate.gif
}


# Combine all data into a single file and plot - useful for when only looking at time based variations
#??ds no use if there are any spatial variations since it plots everything
#??ds have to set titles and columns to be used manually
#??ds fix to use input properly
#??ds should maybe use something other than gif to save?
########################################################################
function plot-time
{
    # Cat the solution files together into a time data file
    cat RESLT/soln* > RESLT/time_soln

    # Set options
    echo "clear" > .plot.gpi
    echo "reset" >> .plot.gpi
    echo "set terminal gif" >> .plot.gpi
    echo "set output \"RESLT/time-plot.gif\"" >> .plot.gpi # ouput file name

    # set ranges of plots
    #echo "set xrange [0:1]" >> .plot.gpi
    #echo "set yrange [-2:2]" >> .plot.gpi


    # plot time data file
    echo "plot \"time_soln\" using 1:2 t \"m\" " >> .plot.gpi


    # Run the commands and leave gnuplot open
    gnuplot -persist .plot.gpi

}