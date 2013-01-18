#!/bin/bash

set -o nounset
set -o errexit

# move to results dir
cd loop_results

# Get end state data
find -name "comp*" | xargs -n1 tail -1 > tmp_data

# Get headers from the first data file combined with some input here:
echo "eps dt mconstraintmethod damping hk" > tmp_filename_headers
find -name "comp*" | head -1 | xargs -n1 head -1 | paste tmp_filename_headers - > tmp_headers

# Convert filenames to input data
find -name "parameters" | xargs -n1 cat > tmp_input_data
paste tmp_input_data tmp_data | \
cat tmp_headers - > ../complete_data

# Get the values at some specified time
time=1
find -name "comp*" | xargs -n1 grep "^$time " | paste tmp_input_data - > ../data_at_time_$time

rm -f tmp_data tmp_input_data tmp_filename_headers tmp_headers
