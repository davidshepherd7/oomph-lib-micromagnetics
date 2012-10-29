#set -o errexit
set -o nounset

make

# Generate the input to tetgen to make a spherical mesh
./generate_tetgen_sphere_input 1.0 1

# Generate mesh, -Y = don't split facets, -q = make a good qualit
    #tetgen -q1.0pa0.01 ./mesh.poly
tetgen -q1.0a0.1p ./mesh.poly
