#! /bin/bash

set -o errexit
set -o nounset

base1=`basename $1`
base2=`basename $2`

echo "format short e" > octave_script.m

echo "load \"$1\";" >> octave_script.m
echo "load \"$2\";" >> octave_script.m

echo "max_diff = max(max(abs( $base1 - $base2)))" >> octave_script.m
echo "mean_diff = mean(mean(abs($base1 - $base2)))" >> octave_script.m

# Run the script
octave -q --persist octave_script.m
