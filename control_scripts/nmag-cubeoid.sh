#!/bin/bash

set -o errexit
set -o nounset

cd semi_implicit_mm_driver

make -k -s

randstring=$(date | md5sum | tr -d ' '| tr -d '-')


# Make sure directory exists and is empty
outdir="../../experiments/nmag_cubeoid_${randstring}"
touch $outdir
mv $outdir ${outdir}.bak
mkdir -p $outdir

# The run command
./semi_implicit_mm_driver -mesh sq_cubeoid -ref 4 -tol 1e-4 -dt 1e-4 -tmax 40.0 -initm xz -happ zero -solver gmres-amg -mag-params nmag-cubeoid -outdir $outdir | tee > $outdir/stdout
