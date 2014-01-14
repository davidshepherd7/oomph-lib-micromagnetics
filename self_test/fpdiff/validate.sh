#! /bin/sh
mkdir -p Validation
./validate.py 2>&1 > Validation/validation.log
