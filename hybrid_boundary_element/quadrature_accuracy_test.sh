#! /bin/bash

set -o errexit
set -o nounset

# The highest Gauss order used by the program
max_i=29

# Number of elements to use in each direction
n_ele=50

# Build and run the program
make
./hybrid_boundary_element_driver $n_ele $n_ele

# Now create the octave script:

echo "format short e" > octave_script.m

# For each possible gauss order:
for (( i=2; i<=$max_i; i++ )); do
    # Load this matrix into octave
    echo "load \"results/boundary_matrix_$i\";" >> octave_script.m
done

# Compare each matrix with the previous (using octave)
for (( i=2; i<=$max_i; i++ )); do
    echo "max_diff($i) = max(max(abs(boundary_matrix_$i - boundary_matrix_$max_i)))" >> octave_script.m
    echo "mean_diff($i) = mean(mean(abs(boundary_matrix_$i - boundary_matrix_$max_i)))" >> octave_script.m
done

# load results from gauss-legendre quadrature
echo "load \"results/gl_50x50_results\"" >> octave_script.m

echo "hold on" >> octave_script.m
echo "semilogy(max_diff(2:end))" >> octave_script.m
echo "semilogy(gl_max_diff(2:end),'r')" >> octave_script.m

echo "figure" >> octave_script.m
echo "hold on" >> octave_script.m
echo "semilogy(mean_diff(2:end))" >> octave_script.m
echo "semilogy(gl_mean_diff(2:end),'r')" >> octave_script.m

# Run the script
octave -q --persist octave_script.m
