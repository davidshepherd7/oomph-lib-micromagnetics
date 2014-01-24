#! /bin/bash

set -o errexit
set -o nounset

mkdir -p Validation

./node_lookup_test

exit 0
