#!/bin/bash

set -o errexit
set -o nounset

make

# Generate the input to tetgen to make a spherical mesh
./generate_tetgen_sphere_input

# Generate mesh, -Y = don't split facets, -q = make a good quality mesh, -p
# = use poly file for input, -a = maximum tet volume
tetgen -q1.0a0.1pY ./mesh.poly

# Run the test driver
./bem_field_test

# Returns failure status if anything fails

exit 0
