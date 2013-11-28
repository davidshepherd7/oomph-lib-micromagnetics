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
cd $CONTROL_SCRIPTS/semi_implicit_mm_driver/
./semi_implicit_mm_driver -dt 0.001 -tmax 0.01 -mesh ut_cubeoid -ref 2 \
    -solver superlu -happ zero -initm xz -outdir $CUBEOID_DIR \
    -mag-params 'simple-llg' \
    > $CUBEOID_DIR/stdout

# Extract + check energies
final_energy $CUBEOID_DIR/trace > $CUBEOID_DIR/energies
wrapped_fpdiff $CUBEOID_DIR/energies $TPWD/../cubeoid_energies
