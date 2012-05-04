set -o errexit
set -o nounset

make

# move old files to backup
touch results/soln0.dat mesh.poly # make sure exist to avoid error in mv
mv -f results/* mesh.* results_old/

# Generate the input to tetgen to make a spherical mesh
./generate_tetgen_sphere_input 1.0 1

# Generate mesh, -Y = don't split facets, -q = make a good qualit
tetgen -q1.0pa0.02 ./mesh.poly

# Run the actual code
./sphere_driver "mesh.1.node" "mesh.1.ele" "mesh.1.face"

# paraview plot
oomph-convert results/soln*.dat
paraview --data="results/soln..vtu" --script="../paraview_vectors.py"
