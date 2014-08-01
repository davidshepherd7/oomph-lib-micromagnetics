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
MM_LIB_DIR="../../"

# Make sure it's all ready
new_clean_dir Validation



# bdf implementation of midpoint
# ============================================================

llg_bdf_midpoint_dir="$TPWD/Validation/llg_bdf_midpoint"
new_clean_dir $llg_bdf_midpoint_dir
cd $CONTROL_SCRIPTS/driver/
./driver llg -ms-method disabled -tmax 0.3 -dt 1e-2 -solver superlu -h-app minus_z \
    -ts midpoint-bdf -newton-tol 1e-12 -initial-m z -mag-params 'simple-llg' \
    -outdir $llg_bdf_midpoint_dir > $llg_bdf_midpoint_dir/stdout

# Check the errors are small by comparing with a file full of zeros
cd $TPWD
cut -d\; -f 4 < $llg_bdf_midpoint_dir/trace > $llg_bdf_midpoint_dir/time_error_norms
wrapped_fpdiff $llg_bdf_midpoint_dir/time_error_norms validata/zeros 0 1e-5

# Standard implementation of midpoint
# ============================================================

llg_midpoint_dir="$TPWD/Validation/llg_midpoint"
new_clean_dir $llg_midpoint_dir
cd $CONTROL_SCRIPTS/driver/
./driver llg -ms-method disabled -tmax 0.3 -dt 1e-2 -solver superlu -h-app minus_z \
    -ts old-imr -newton-tol 1e-12 -initial-m z -mag-params 'simple-llg' \
    -outdir $llg_midpoint_dir > $llg_midpoint_dir/stdout

# Check the errors are small by comparing with a file full of zeros
cd $TPWD
cut -d\; -f 4 < $llg_midpoint_dir/trace > $llg_midpoint_dir/time_error_norms
wrapped_fpdiff $llg_midpoint_dir/time_error_norms validata/zeros 0 1e-5


# same for ll residual
# ============================================================


# Run simulation with ll residual
ll_dir="$TPWD/Validation/ll"
new_clean_dir $ll_dir
cd $CONTROL_SCRIPTS/driver/
./driver llg -ms-method disabled -dt 1e-2 -tmax 0.3 -solver superlu -h-app minus_z \
    -ts imr -initial-m z -mag-params 'simple-llg' -fd-jac \
    -outdir $ll_dir > $ll_dir/stdout

# Check the errors are small
cd $TPWD
cut -d\; -f 4 < $ll_dir/trace > $ll_dir/time_error_norms
wrapped_fpdiff $ll_dir/time_error_norms validata/zeros 0 1e-5



# llg residual with two separate squares of mesh
# ============================================================

# Run simulation with llg residual. Midpoint is ESSENTIAL to get low enough
# errors without tiny step size
mul_mesh_dir="$TPWD/Validation/mul_mesh_llg"
new_clean_dir $mul_mesh_dir
cd $CONTROL_SCRIPTS/driver/
./driver llg -ms-method disabled -tmax 0.3 -dt 1e-2 -solver superlu -h-app minus_z \
    -ts imr -newton-tol 1e-12 -mesh multi_ut_square \
    -initial-m z -mag-params 'simple-llg' \
    -outdir $mul_mesh_dir > $mul_mesh_dir/stdout

# Check the errors are small by comparing with a file full of zeros
cd $TPWD
cut -d\; -f 4 < $mul_mesh_dir/trace > $mul_mesh_dir/time_error_norms
wrapped_fpdiff $mul_mesh_dir/time_error_norms validata/zeros 0 1e-5


# Explicit timestepping
# ============================================================
explicit_dir="$TPWD/Validation/explicit_llg"
new_clean_dir $explicit_dir
cd $CONTROL_SCRIPTS/driver/
./driver ll -ms-method disabled -tmax 0.3 -dt 1e-2 -h-app minus_z \
    -ts rk4 -mesh q_single_element \
    -initial-m z -mag-params 'simple-llg' \
    -outdir $explicit_dir > $explicit_dir/stdout

# Check the errors are small by comparing with a file full of zeros
cd $TPWD
cut -d\; -f 4 < $explicit_dir/trace > $explicit_dir/time_error_norms
wrapped_fpdiff $explicit_dir/time_error_norms validata/zeros 0 1e-5


exit 0
