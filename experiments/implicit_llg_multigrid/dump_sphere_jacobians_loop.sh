#!/bin/bash

set -o nounset

program_name=`pwd`"/multigrid_experiments"
BASE=`pwd`"/runs/sphere_sin_cos_initial"

refines="1 2 3"
dts="0.1 0.01 0.001"

# refines="1"
# dts="0.1"

make

for refine in $refines; do
    for dt in $dts; do
        DIR="${BASE}/refine${refine}_dt${dt}"
        echo $DIR
        mkdir -p $DIR
        ./multigrid_experiments '-dt' $dt '-tmax' $dt '-outdir' $DIR > ${DIR}/trace
        mv "${DIR}/jacobian_t0" "${DIR}/jacobian_t0_refine${refine}_dt${dt}"
    done
done