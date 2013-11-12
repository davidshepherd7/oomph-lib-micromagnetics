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

new_clean_dir $TPWD/Validation
