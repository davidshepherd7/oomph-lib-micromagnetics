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


# llg residual
# ============================================================

# Run simulation with llg residual. Midpoint is ESSENTIAL to get low enough
# errors without tiny step size
llg_dir="$TPWD/Validation/llg"
new_clean_dir $llg_dir
cd $CONTROL_SCRIPTS/llg_driver/
./llg_driver -resi llg -tmax 0.3 -dt 1e-2 -solver superlu -happ minus_z \
    -ts midpoint -newton-tol 1e-12 \
    -initm z -outdir $llg_dir -mag-params 'simple-llg' \
    > $llg_dir/stdout

# Check the errors are small by comparing with a file full of zeros
cd $TPWD
cut -d\; -f 4 < $llg_dir/trace > $llg_dir/time_error_norms
wrapped_fpdiff $llg_dir/time_error_norms validata/zeros 0 1e-5


# same for ll residual
# ============================================================


# Run simulation with ll residual
ll_dir="$TPWD/Validation/ll"
new_clean_dir $ll_dir
cd $CONTROL_SCRIPTS/llg_driver/
./llg_driver -resi ll -dt 1e-2 -tmax 0.3 -solver superlu -happ minus_z \
    -ts midpoint \
    -initm z -outdir $ll_dir -mag-params 'simple-llg' -fd-jac \
    > $ll_dir/stdout

# Check the errors are small
cd $TPWD
cut -d\; -f 4 < $ll_dir/trace > $ll_dir/time_error_norms
wrapped_fpdiff $ll_dir/time_error_norms validata/zeros 0 1e-5



# llg residual with two separate squares of mesh
# ============================================================

# Run simulation with llg residual. Midpoint is ESSENTIAL to get low enough
# errors without tiny step size
mul_mesh_dir="$TPWD/Validation/llg"
new_clean_dir $mul_mesh_dir
cd $CONTROL_SCRIPTS/llg_driver/
./llg_driver -resi llg -tmax 0.3 -dt 1e-2 -solver superlu -happ minus_z \
    -ts midpoint -newton-tol 1e-12 -mesh multi_ut_square \
    -initm z -outdir $mul_mesh_dir -mag-params 'simple-llg' \
    > $mul_mesh_dir/stdout

# Check the errors are small by comparing with a file full of zeros
cd $TPWD
cut -d\; -f 4 < $mul_mesh_dir/trace > $mul_mesh_dir/time_error_norms
wrapped_fpdiff $mul_mesh_dir/time_error_norms validata/zeros 0 1e-5

exit 0
