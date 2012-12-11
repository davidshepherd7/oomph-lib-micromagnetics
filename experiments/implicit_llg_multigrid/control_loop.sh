#!/bin/bash

set -o nounset
set -o errexit

make

# Initialise output files
echo "nx,dt,smoother,vc,avg_newton_its,avg_solver_its" > iteration_counts


# Run lots of different grid sizes in parallel

./pmg_loop.sh 10 &

./pmg_loop.sh 20 &

./pmg_loop.sh 50 &

./pmg_loop.sh 100 &
