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
MM_LIB_DIR="../../"

# Make sure it's all ready
new_clean_dir Validation

# Run
./mallinson_test > Validation/mallinson_dynamics

# Check
wrapped_fpdiff validata/mallinson_dynamics Validation/mallinson_dynamics

exit 0
