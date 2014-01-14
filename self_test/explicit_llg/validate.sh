#! /bin/sh
mkdir -p Validation
# Run the real script and exit with it's error status
./validate.py 2>&1 > Validation/validation.log
