#!/bin/bash

set -o nounset

make

nxs="5 8 10 12"
nrefines="1 2 3"

# generate meshes
tetgen -q1.2 -p -a400 cubeoid # generates cubeoid.1....
tetgen -q1.2 -r -a200 cubeoid.1
tetgen -q1.2 -r -a100 cubeoid.2
tetgen -q1.2 -r -a50 cubeoid.3
tetgen -q1.2 -r -a25 cubeoid.4  # generates cubeoid.5.....

for nx in $nxs; do
    outdir="struct_${nx}"
    echo $outdir
    rm -r $outdir && mkdir $outdir
    echo "time mx my mz" > $outdir/averages
    echo "time hx hy hz" > $outdir/field_averages
    ./semi_implicit_driver -structured_mesh $nx -outdir $outdir > "${outdir}/trace"
done &

for nrefine in $nrefines; do
    outdir="unstruct_$nrefine"
    rm -r $outdir && mkdir $outdir
    echo $outdir
    echo "time mx my mz" > $outdir/averages
    echo "time hx hy hz" > $outdir/field_averages
    ./semi_implicit_driver -unstructured_mesh $nrefine -outdir $outdir > "${outdir}/trace"
done &
