#!/bin/bash

set -o errexit
set -o nounset

make -k

# as="3 4 5" # 0.001 0.0001"

# # rebuild mesh files
# #fig2poly circle.fig
# triangle -q20 circle.fig.poly
# triangle -rp -q30 -a1.6 circle.fig.1.poly
# triangle -rp -q30 -a0.4 circle.fig.2.poly
# triangle -rp -q30 -a0.1 circle.fig.3.poly
# triangle -rp -q30 -a0.025 circle.fig.4.poly
# triangle -rp -q30 -a0.00625 circle.fig.5.poly


# # rebuild corner file
# ./poly2corner.sh circle.fig.poly

# for a in $as; do

#     ./circle_bem -mesh "circle.fig.${a}" -corners "circle.fig.poly.corner"
    
#     oomph-convert results/*.dat

#     mkdir -p results_${a}
#     mv results/* results_${a}
    
# done
 
rm -rf results/*

./circle_bem -mesh "circle.fig.3" -corners "circle.fig.poly.corner"

cat results/errors

 #oomph-convert results/soln*.dat

#paraview --data="results/soln..vtu" &
