#!/bin/sh

set -o errexit

# Run the python script and return its error status
./check_jacobians.py
