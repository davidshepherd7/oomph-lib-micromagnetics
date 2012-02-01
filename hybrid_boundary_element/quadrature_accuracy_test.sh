#! /bin/bash

set -o errexit
set -o nounset

# For a few different numbers of elements
for j in 5 10 50 250; do

    # Create a script for octave to parse the data
    echo "disp($j)" > octave_script.m

    # For each possible gauss order:
    for i in 2 3 4; do
    # Horrible hack to change gauss order in script:
	echo "const unsigned gauss_order = $i;" > gauss_order.h

    # Build and run the program
	make > /dev/null
	./hybrid_boundary_element_driver $j $j > /dev/null

    # Move the resulting file
	mv "results/boundary_matrix" "results/boundary_matrix_$i"

    # Load this matrix into octave
	echo "load \"results/boundary_matrix_$i\";" >> octave_script.m
    done

    # Compare each matrix with the previous (using octave)
    for i in 3 4; do
	echo "max(max(boundary_matrix_$i - boundary_matrix_$(($i-1))))" >> octave_script.m
    done

    octave -q octave_script.m

done
