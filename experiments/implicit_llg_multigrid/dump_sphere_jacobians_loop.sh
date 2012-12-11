#!/bin/bash

set -o nounset

program_name=`pwd`"/multigrid_experiments"
BASE=`pwd`"/runs/sphere"

refines="1 2 3"
dts="0.1 0.01 0.001"

# generate tetgen meshes (do by hand since bash has no decimal maths)
tetgen sphere.poly -qa0.005
tetgen -qra0.001 sphere.1
tetgen -qra0.00025 sphere.2
tetgen -qra0.00005 sphere.3
#tetgen -qra0.000001 sphere.4

# refines="1"
# dts="0.1"

make

for refine in $refines; do
    for dt in $dts; do
        DIR="${BASE}/refine${refine}_dt${dt}"
        echo $DIR
        mkdir -p $DIR
        ./multigrid_experiments '-dt' $dt '-tmax' 0 '-outdir' $DIR '-r' $refine '-mesh' 'sphere' > ${DIR}/trace
        #mv "${DIR}/jacobian_t0" "${DIR}/jacobian_t0_refine${refine}_dt${dt}"
    done
done