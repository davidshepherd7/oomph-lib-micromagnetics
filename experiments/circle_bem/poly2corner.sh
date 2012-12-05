#!/bin/bash

set -o errexit
set -o nounset

# Extract x and y coords for nodes from a .poly file

# remove first line, remove comments, remove all after first blank line (i.e. only care about nodes) then take only fields 2 and 3 from what remains (i.e. x and y coords).
cat $1 | sed '1d' | sed 's/#.*//g' | sed -n '/^$/q;p' | cut -f2,3 -d' ' | perl -pe 'chomp if eof' > $1.corner