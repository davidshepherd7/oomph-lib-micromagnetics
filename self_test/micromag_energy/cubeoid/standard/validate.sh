#!/bin/sh

micromagroot="../../../.."
CONTROL_SCRIPTS="$micromagroot/control_scripts"
MM_LIB_DIR="$micromagroot/"

# source common commands for energy tests
. "$micromagroot/self_test/micromag_energy/init_validate.sh"


# Cubeoid test
# ============================================================

# This gives us some finite exchange and magnetostatic energies to check
# against previous runs.

CUBEOID_DIR=$TPWD/Validation
new_clean_dir $CUBEOID_DIR

# Run simulation
cd $CONTROL_SCRIPTS/driver llg -ms-method decoupled/
./driver llg -ms-method decoupled -dt 0.001 -tmax 0.01 -mesh ut_cubeoid -ref 3 \
    -solver superlu -h-app zero -initial-m xz -mag-params 'simple-llg' \
    -hlib-bem 0 \
    -outdir $CUBEOID_DIR > $CUBEOID_DIR/stdout

# Extract + check energies
final_energy $CUBEOID_DIR/trace > $CUBEOID_DIR/energies
wrapped_fpdiff $CUBEOID_DIR/energies $TPWD/../cubeoid_energies 8 1.5e-4
