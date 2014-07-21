#!/bin/bash

set -o errexit
set -o nounset


# convert code to something that will compile under c++98, should compile
# but might make comments not make sense...

omm_root="$(pwd)/../"

code_files=$(find "$omm_root" -name '*.h' -o -name '*.cc')

echo "c++ files found: $code_files"

sed -e 's/ override//g' -e 's/unique_ptr/auto_ptr/g' -i $code_files
