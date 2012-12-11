#!/bin/bash

set -o nounset

# function parallel_exec()
# {
#     while $jobsrunning < $maxjobs
#     {
#         do $1 &
#         $jobsrunning += 1
#     }
#     wait
# }

function loop_core()
{
    
    echo $LABEL
    mkdir -p $DIR
    ./multigrid_experiments '-nx' $nx '-ny' $nx '-dt' $dt '-tmax' $dt \
        '-amgsmth' $smoother '-amgvcyc' $vc '-outdir' $DIR \
        > ${DIR}/trace
    mv "${DIR}/jacobian_t0" "${DIR}/jacobian_${LABEL}"
}

program_name=`pwd`"/multigrid_experiments"
BASE=`pwd`"/runs/multigrid"

rm -rf $BASE/*

nxs="20 50 100"
dts="1e-1 1e-2 1e-4 1e-6"

# 0 = Jacobi
# 1 = Gauss-Seidel, sequential
#     (very slow in parallel!)
# 2 = Gauss-Seidel, interior points in parallel, boundary sequential
#     (slow in parallel!)
# 3 = hybrid Gauss-Seidel or SOR, forward solve
# 4 = hybrid Gauss-Seidel or SOR, backward solve
# !!6 = hybrid symmetric Gauss-Seidel or SSOR - overlaps with complex so disabled for simplicity
# 6 = Schwarz
# 7 = Pilut
# 8 = ParaSails
# 9 = Euclid

smoothers="0 1 2 3 4 5 6 7 8 9"

v_cycles="1 2"

make
for smoother in $smoothers; do

    # iteration_counts file header
    echo > iteration_counts_$smoother
    echo "# nxs = $nxs" >> iteration_counts_$smoother
    echo "# dts = $dts" >> iteration_counts_$smoother
    echo "# Missing values indicate failure to finish the run." >> iteration_counts_$smoother
    echo "nx dt vc avg_newton_its avg_solver_its" >> iteration_counts_$smoother

    for vc in $v_cycles; do
        for dt in $dts; do 
            cat nxs | parallel 'DIR="${BASE}/nx${nx}_dt${dt}_smth${smoother}_vc_$vc"; mkdir -p $DIR; nx={}; ./multigrid_experiments "-nx" $nx "-ny" $nx "-dt" $dt "-tmax" $dt "-amgsmth" $smoother "-amgvcyc" $vc "-outdir" $DIR > ${DIR}/trace; mv "${DIR}/jacobian_t0" "${DIR}/jacobian_${LABEL}"'
        done
    done

    # print the table of results
    cat iteration_counts_$smoother | columns -s, -t

done

