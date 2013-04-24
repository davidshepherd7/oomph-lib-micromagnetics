#!/bin/bash

set -o errexit
set -o nounset

make

./llg_exchange_test

exit fpdiff.py "./validation/precomputed_solution_at_t0.03.dat" "./validation/generated_solution_at_t0.03.dat"
