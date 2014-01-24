#!/bin/bash

set -o errexit
set -o nounset

mkdir -p Validation

./sum_of_matrices_test > Validation/sum_matrices_test_trace

exit 0
