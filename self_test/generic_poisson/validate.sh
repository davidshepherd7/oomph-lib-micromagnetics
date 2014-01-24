#!/bin/bash

set -o errexit
set -o nounset

mkdir -p Validation

./generic_poisson_test
exit 0
