#set -o errexit
set -o nounset

make

sort_clean_file()
{
    for file in $@; do
	sort $file -k 1,1n -k 2,2n | grep -v -E -e " 0$" -e "e-[1-3][0-9]" -e "e-0[7-9]" > ${file}_sorted

	mv ${file}_sorted $file
	rm -f ${file}_sorted
    done
}


# move old files to backup
touch results/soln0.dat mesh.poly # make sure exist to avoid error in mv
mv -f results/* mesh.* results_old/
rm -rf matrices/
mkdir matrices

rm -rf fd non-fd
mkdir fd non-fd

# Generate the input to tetgen to make a spherical mesh
./generate_tetgen_sphere_input 1.0 1

# Generate mesh, -Y = don't split facets, -q = make a good qualit
    #tetgen -q1.0pa0.01 ./mesh.poly
tetgen -q1.0a0.1p ./mesh.poly

#tetgen ./mesh.poly -qp

# Run the actual code
./sphere_driver "mesh.1.node" "mesh.1.ele" "mesh.1.face"

# # paraview plot
# oomph-convert results/soln*.dat && paraview --data="results/soln..vtu"


# Clean up all output
sort_clean_file ./fd/* ./non-fd/*

cd "./non-fd/"
for file in *; do
    if [ -f "../fd/$file" ] ; then
	echo -e "Taking diff of $file\n"
	fpdiff.py "$file" "../fd/$file" 0.1 1e-7 | tee "../${file}_fp_diff"
    else
	echo "No finite differenced equivalent of $file"
    fi
done
