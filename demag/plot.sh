#!/bin/bash

# use -persist to keep plot window open permenantly

# Set options
echo "clear" > .plot.gpi
echo "reset" >> .plot.gpi
echo "set terminal gif animate delay 20" >> .plot.gpi
echo "set output \"RESLT/animate.gif\"" >> .plot.gpi
echo "set isosample 40" >> .plot.gpi
echo "set hidden3d" >> .plot.gpi
echo "set xrange [0:1]" >> .plot.gpi
echo "set yrange [-1:1]" >> .plot.gpi


# plot data files
for file in $@
do
    echo "plot \"$file\" using 1:3, \"$file\" using 1:4, \"$file\" using 1:5" >> .plot.gpi
done

# Run the commands
gnuplot .plot.gpi

# Clean up
#rm .plot.gpi

# Show the gif
eog RESLT/animate.gif