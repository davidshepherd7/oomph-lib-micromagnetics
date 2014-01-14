#! /bin/sh
set -o errexit
mkdir -p "Validation"
./hlib_matrix_generate_driver 2>&1 > "Validation/validation.log"

fpdiff.py "Validation/new_bem_matrix" "Validation/old_bem_matrix" 1 1e-10
