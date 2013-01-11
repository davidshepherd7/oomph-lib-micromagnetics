#!/bin/bash

set -o nounset
set -o errexit

# Get end state data
find -name "comp*" | xargs -n1 tail -1 > tmp_data

# Get headers
echo "eps renormalise damping hk" > tmp_filename_headers
find -name "comp*" | head -1 | xargs -n1 head -1 | paste tmp_filename_headers - > tmp_headers

# Convert filenames to input data
find -name "comp*" | xargs -n1 dirname | sed 's|./loop_results/sp_||' | \
tr _ ' ' | tr 'm' '-' | \
paste - tmp_data | \
cat tmp_headers - > complete_data

rm -f tmp_data  tmp_filename_headers tmp_headers
