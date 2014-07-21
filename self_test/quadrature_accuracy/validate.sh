
#!/bin/sh
set -o errexit
set -o nounset

# Clean Validation dir
touch Validation
rm -r Validation
mkdir Validation

# Run the binary
./quadrature_accuracy_driver
