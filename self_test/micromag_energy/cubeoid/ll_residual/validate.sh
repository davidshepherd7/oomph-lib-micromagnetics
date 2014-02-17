#!/bin/sh

micromagroot="../../../.."
CONTROL_SCRIPTS="$micromagroot/control_scripts"
MM_LIB_DIR="$micromagroot"

# source common commands for energy tests
. "$micromagroot/self_test/micromag_energy/init_validate.sh"


# Landau-Lifshitz residual cubeoid test
# ============================================================

# This gives us some finite exchange and magnetostatic energies to check
# against previous runs. Should be the same as for llg.

LL_CUBEOID_DIR=$TPWD/Validation
new_clean_dir $LL_CUBEOID_DIR

# Run simulation
cd $CONTROL_SCRIPTS/driver llg -ms-method decoupled/
./driver ll -ms-method decoupled -dt 0.001 -tmax 0.01 -mesh ut_cubeoid -ref 3 \
    -solver superlu -h-app zero -initial-m xz -mag-params 'simple-llg' \
    -fd-jac -hlib-bem 0 \
    -outdir $LL_CUBEOID_DIR> $LL_CUBEOID_DIR/stdout

# Extract + check energies
final_energy $LL_CUBEOID_DIR/trace > $LL_CUBEOID_DIR/energies
wrapped_fpdiff $LL_CUBEOID_DIR/energies $TPWD/../cubeoid_energies 8 1.5e-4
