#!/bin/sh

micromagroot="../../../.."
CONTROL_SCRIPTS="$micromagroot/control_scripts"
MM_LIB_DIR="$micromagroot"

# source common commands for energy tests
. "$micromagroot/self_test/micromag_energy/init_validate.sh"

# Two close spheres test
# ============================================================

# Now test that the energy is ~ doubled when we run for two well separated
# spheres.

TWO_SPHERE_DIR=$TPWD/Validation
new_clean_dir $TWO_SPHERE_DIR

# Run simulation
cd $CONTROL_SCRIPTS/semi_implicit_mm_driver/
./semi_implicit_mm_driver -dt 0.001 -tmax 0.001 -ref 3 \
    -solver superlu -happ minus_z -initm exactly_z \
    -mag-params 'simple-llg' -mesh multi_ut_sphere -xshift 1.1 \
    -outdir $TWO_SPHERE_DIR > $TWO_SPHERE_DIR/stdout

# Extract + check energies
final_energy $TWO_SPHERE_DIR/trace > $TWO_SPHERE_DIR/energies
wrapped_fpdiff $TWO_SPHERE_DIR/energies $TPWD/../two_sphere_energies 15 1e-4

# ??ds this should fail for now
