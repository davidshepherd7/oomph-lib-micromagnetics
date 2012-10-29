#!/bin/bash

set -o nounset

program_name=`pwd`"/multigrid_experiments"
BASE=`pwd`"/runs/sin_cos_initial"

nxs="10 20 40"
dts="0.1 0.01 0.001"

# nxs="10"
# dts="0.1"

make

for nx in $nxs; do
    for dt in $dts; do
        DIR="${BASE}/nx${nx}_dt${dt}"
        echo $DIR
        mkdir -p $DIR
        ./multigrid_experiments '-nx' $nx '-ny' $nx '-dt' $dt '-tmax' $dt '-outdir' $DIR > ${DIR}/trace
        mv "${DIR}/jacobian_t0" "${DIR}/jacobian_t0_nx${nx}_dt${dt}"
    done
done