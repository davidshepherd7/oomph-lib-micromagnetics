#!/bin/bash

set -o errexit

PATH="$PATH:$HOME/oomph-lib/bin"

make

rm -r results/*

./multigrid_experiments

# pretty print last 10 results with headings
cat <(head iteration_counts -n1)  <(tail iteration_counts) | column -s, -t

# oomph-convert results/*.dat
# paraview results/soln..vtu &
