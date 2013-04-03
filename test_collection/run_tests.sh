#!/bin/bash

set -o errexit

# run tests using oomph_root/bin/parallel_self_test.py, eg.
../../../bin/parallel_self_test.py -C ../ $@
