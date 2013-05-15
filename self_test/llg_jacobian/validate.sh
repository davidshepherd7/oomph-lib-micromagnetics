#!/bin/sh

# Run the python script and return its error status
./check_jacobians.py --fast --serial
exit $?
