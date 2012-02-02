#! /bin/bash

set -o errexit
set -o nounset

rm -f octave_script.m

# # For a few different numbers of elements
## for j in 5 10 50 250; do

max_i=29

    # Build and run the program
make
./hybrid_boundary_element_driver 5 5

    # For each possible gauss order:
for (( i=2; i<=$max_i; i++ )); do
    # Load this matrix into octave
    echo "load \"results/boundary_matrix_$i\";" >> octave_script.m
done

    # Compare each matrix with the previous (using octave)
for (( i=2; i<=$max_i; i++ )); do
    echo "results($i) = max(max(abs(boundary_matrix_$i - boundary_matrix_$max_i)));" >> octave_script.m
done

octave --persist octave_script.m

# done
