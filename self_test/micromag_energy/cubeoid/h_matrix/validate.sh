#!/bin/sh

micromagroot="../../../.."
CONTROL_SCRIPTS="$micromagroot/control_scripts"
MM_LIB_DIR="$micromagroot"

# source common commands for energy tests
. "$micromagroot/self_test/micromag_energy/init_validate.sh"

# This gives us some finite exchange and magnetostatic energies to check
# against previous runs. Should be the same as for llg.

CUBEOID_DIR=$TPWD/Validation
new_clean_dir $CUBEOID_DIR

# Run simulation
cd $CONTROL_SCRIPTS/driver/
./driver llg -decoupled-ms -dt 0.001 -tmax 0.001 -mesh ut_cubeoid  -ref 2 \
    -solver superlu -happ zero -initm xz -mag-params 'simple-llg' \
    -hierarchical_bem 1 \
     -outdir $CUBEOID_DIR 2>&1 > $CUBEOID_DIR/stdout

# Extract + check energies
final_energy $CUBEOID_DIR/trace > $CUBEOID_DIR/energies
wrapped_fpdiff $CUBEOID_DIR/energies $TPWD/../cubeoid_energies 15 1e-4
