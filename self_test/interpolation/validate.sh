#!/bin/bash

set -o errexit
set -o nounset

mkdir -p "Validation"

./polynomial_interpolation_test

exit 0
