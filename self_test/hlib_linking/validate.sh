#! /bin/sh
set -o errexit
mkdir -p "Validation"
./hlib_test_driver 2>&1 > "Validation/validation.log"
