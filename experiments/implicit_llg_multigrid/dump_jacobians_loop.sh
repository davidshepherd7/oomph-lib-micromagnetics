#!/bin/bash

set -o nounset

program_name=`pwd`"/multigrid_experiments"
SBASE=`pwd`"/runs/structured_squares"

nxs="20 40 80 160"
dts="0.1 0.05 0.01 0.001"
field_dirs="x y z"

# nxs="10"
# dts="0.1"

make

for nx in $nxs; do
    for dt in $dts; do
        for field_dir in $field_dirs; do
            DIR="${SBASE}/nx${nx}_dt${dt}_f${field_dir}"
            echo $DIR
            mkdir -p $DIR

            ./multigrid_experiments '-nx' $nx '-ny' $nx '-dt' $dt \
                '-tmax' $dt  '-outdir' $DIR '-field' $field_dir > ${DIR}/trace
            mv "${DIR}/jacobian_tmax" "${DIR}/jacobian_tmax_nx${nx}_dt${dt}_f${field_dir}"
        done
    done
done


# Now do unstrctured grids
# ============================================================

UBASE=`pwd`"/runs/unstructured_squares"
refines="1 2 3 4 5"

triangle -q20 -a0.1 square.poly
triangle -rp -q30 -a0.025 square.1.poly
triangle -rp -q30 -a0.01 square.2.poly
triangle -rp -q30 -a0.0025 square.3.poly
triangle -rp -q30 -a0.0001 square.4.poly
triangle -rp -q30 -a0.000025 square.5.poly

for refine in $refines; do
    for dt in $dts; do
        for field_dir in $field_dirs; do
            DIR="${UBASE}/refine${refine}_dt${dt}_f${field_dir}"
            echo $DIR
            mkdir -p $DIR
            ./multigrid_experiments '-r' $refine '-dt' $dt '-tmax' $dt \
                '-outdir' $DIR '-field' $field_dir > ${DIR}/trace
            mv "${DIR}/jacobian_tmax" "${DIR}/jacobian_tmax_refine${refine}_dt${dt}_f${field_dir}"
        done
    done
done
