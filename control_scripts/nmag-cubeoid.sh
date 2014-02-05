#!/bin/bash

set -o errexit
set -o nounset

make -k -s -C ./semi_implicit_mm_driver/

randstring=$(date | md5sum | tr -d ' '| tr -d '-')


# Make sure directory exists and is empty
outdir="../experiments/nmag_cubeoid_${randstring}"
mkdir -p $outdir

echo "Output going into $outdir"

# The run command
./semi_implicit_mm_driver/semi_implicit_mm_driver -mesh sq_cubeoid -ref 3 -tol 1e-4 -dt 1e-4 -tmax 40.0 -initm xz -happ zero -solver gmres -preconditioner amg -mag-params nmag-cubeoid -outdir $outdir 2>&1 | tee $outdir/stdout
