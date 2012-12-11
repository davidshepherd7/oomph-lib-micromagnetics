#!/bin/bash

set -o nounset

make

nxs="5 8 10 12"

for nx in $nxs; do
    outdir="quad_${nx}"
    echo $outdir
    rm -rf $outdir && mkdir $outdir
    echo "time mx my mz" > $outdir/averages
    echo "time hx hy hz" > $outdir/field_averages
    ./semi_implicit_driver -nx $nx -outdir $outdir > "${outdir}/trace"
done &
