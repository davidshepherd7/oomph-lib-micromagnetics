#!/bin/bash

set -o errexit
set -o nounset

make

touch results && rm -rf results &&  mkdir results

echo "time mx my mz" > results/averages
echo "time hx hy hz" > results/field_averages

# tetgen -qa30 cubeoid.poly

./semi_implicit_driver
