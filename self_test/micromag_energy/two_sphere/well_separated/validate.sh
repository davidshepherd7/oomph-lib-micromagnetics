#!/bin/sh

micromagroot="../../../.."
CONTROL_SCRIPTS="$micromagroot/control_scripts"
MM_LIB_DIR="$micromagroot"

# source common commands for energy tests
. "$micromagroot/self_test/micromag_energy/init_validate.sh"


# Two separate sphere test
# ============================================================

# Now test that the energy is ~ doubled when we run for two well separated
# spheres.

TWO_SPHERE_DIR=$TPWD/Validation
new_clean_dir $TWO_SPHERE_DIR

# Run simulation
cd $CONTROL_SCRIPTS/driver llg -decoupled-ms/
./driver llg -decoupled-ms -dt 0.001 -tmax 0.001 -ref 3 \
    -solver superlu -h-app minus_z -initial-m exactly_z  \
    -mag-params 'simple-llg' -mesh multi_ut_sphere -xshift 20 \
    -hierarchical-bem 0 \
    -outdir $TWO_SPHERE_DIR > $TWO_SPHERE_DIR/stdout

# Extract + check energies
final_energy $TWO_SPHERE_DIR/trace > $TWO_SPHERE_DIR/energies
wrapped_fpdiff $TWO_SPHERE_DIR/energies $TPWD/../two_sphere_energies 15 1e-4



# If M_s=1, r=1, H_app=[0,0,-1.1] (numbers from copying formula into python)
# E_ms = -0.5 * 1 * (-0.33333 * 1) * (4.0/3) * pi = 0.6981247194807239
# E_zee = -1 * -1.1 * (4.0/3) * pi = 4.60766922526503
# E_ca = 0 (because the magnetics is aligned with easy axis)


# Two well separated spheres have double this energy

# Unfortunately good approximations to a sphere are EXPENSIVE.... so it's
# only roughly right :( Hence 1e-4 as numerical zero and 15% error allowed.
