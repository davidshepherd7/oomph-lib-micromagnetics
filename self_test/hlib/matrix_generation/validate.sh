#! /bin/sh
set -o errexit
mkdir -p "Validation"

# H-bem done outside of problem
./hlib_matrix_generate_driver -mesh ut_cubeoid -ref 3 -outdir Validation \
    2>&1 > "Validation/validation.log"

# ??ds loose tols for now, tighten up later? Or maybe measure eigs instead?
fpdiff.py "Validation/new_bem_matrix" "Validation/old_bem_matrix" 5 5e-5 \
    2>&1 > "Validation/validation.log"
