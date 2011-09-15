#!/bin/bash

# use -persist to keep plot window open permenantly

# Set options
echo "clear" > .plot.gpi
echo "reset" >> .plot.gpi
echo "set terminal gif animate delay 8" >> .plot.gpi # delay = time between frames
echo "set output \"RESLT/animate.gif\"" >> .plot.gpi # ouput file name

# set ranges of plots
echo "set xrange [0:1]" >> .plot.gpi
echo "set yrange [-1:1]" >> .plot.gpi


# plot data files
for file in $@
do
    echo "plot \"$file\" using 1:3 title \"M_x\", \"$file\" using 1:4 t \"M_y\", \"$file\" using 1:5 t \"M_z\" " >> .plot.gpi
done

# Run the commands
gnuplot .plot.gpi

# Show the gif
eog RESLT/animate.gif