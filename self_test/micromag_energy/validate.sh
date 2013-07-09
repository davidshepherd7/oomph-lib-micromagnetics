#!/bin/sh

set -o errexit

HAVE_PYTHON=$1

set -o nounset

new_clean_dir()
{
    touch $1
    rm -r $1
    mkdir $1
}

final_energy()
{
    (head -1 < $1; tail -1 < $1) | cut -d';' -f28-31 | tr -d ';'
}

# An fpdiff wrapper that always puts the output in the right place and
# checks if we can run fpdiff.
VALDIR=$(pwd)/Validation
wrapped_fpdiff()
{
    mkdir -p $VALDIR

    if test "$HAVE_PYTHON" = "no_fpdiff"; then
        echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
    else
        ../../../../bin/fpdiff.py $@ >> $VALDIR/validation.log
    fi
}

TPWD=$(pwd)
CONTROL_SCRIPTS="../../control_scripts"

make -k -s -C $CONTROL_SCRIPTS/semi_implicit_mm_driver/
new_clean_dir $TPWD/Validation


# Sphere test
# ============================================================

# This gives us exactly calcuable zeeman, magnetostatic and anisotropy
# energies and zero exchange.

SPHERE_DIR=$TPWD/Validation/sphere
new_clean_dir $SPHERE_DIR

# Run simulation
cd $CONTROL_SCRIPTS/semi_implicit_mm_driver/
./semi_implicit_mm_driver -dt 0.001 -tmax 0.001 -ref 3 -mesh ut_sphere -solver superlu -happ minus_z -initm exactly_z -outdir $SPHERE_DIR -mag-params 'simple-llg' \
    > $SPHERE_DIR/stdout

# Extract + check energies
final_energy $SPHERE_DIR/trace > $SPHERE_DIR/energies
wrapped_fpdiff $SPHERE_DIR/energies $TPWD/validata/sphere_energies 15 1e-4


# If M_s=1, r=1, H_app=[0,0,-1.1] (numbers from copying formula into python)
# E_ms = -0.5 * 1 * (-0.33333 * 1) * (4.0/3) * pi = 0.6981247194807239
# E_zee = -1 * -1.1 * (4.0/3) * pi = 4.60766922526503
# E_ca = 0 (because the magnetics is aligned with easy axis)

# Unfortunately good approximations to a sphere are EXPENSIVE.... so it's
# only roughly right :( Hence 1e-4 as numerical zero and 15% error allowed.



# Cubeoid test
# ============================================================

# This gives us some finite exchange and magnetostatic energies to check
# against previous runs.

CUBEOID_DIR=$TPWD/Validation/cubeoid
new_clean_dir $CUBEOID_DIR

# Run simulation
cd $CONTROL_SCRIPTS/semi_implicit_mm_driver/
./semi_implicit_mm_driver -dt 0.001 -tmax 0.01 -ref 2 -mesh sq_cubeoid -solver superlu -happ zero -initm xz -outdir $CUBEOID_DIR -mag-params 'simple-llg' \
    > $CUBEOID_DIR/stdout

# Extract + check energies
final_energy $CUBEOID_DIR/trace > $CUBEOID_DIR/energies
wrapped_fpdiff $CUBEOID_DIR/energies $TPWD/validata/cubeoid_energies
