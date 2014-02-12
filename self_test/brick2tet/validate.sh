
#!/bin/sh
set -o errexit
set -o nounset

# Clean Validation dir
touch Validation
rm -r Validation
mkdir Validation

mkdir Validation/brick
mkdir Validation/tet

# Run the binary
./brick2tet_driver

# Compare solution at nodes
fpdiff.py "Validation/brick.dump" "Validation/tet.dump" 0.1 1e-8
