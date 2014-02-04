#!/bin/sh

set -o errexit
set -o nounset

touch Validation
rm -r Validation
mkdir Validation

# Run the binary and return its error status
./square2tri_driver

echo "" > Validation/validation.log
fpdiff.py Validation/real_tri_mesh Validation/tri_mesh 2>&1 # >> Validation/validation.log

exit 0
