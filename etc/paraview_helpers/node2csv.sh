#!/bin/bash

cat <(echo "n,x,y,z?") \
    <(sed < $1 -e 's/^ *//g' -e '2d' -e '$d' -e 's/ \+/,/g') \
    | \
    cut -d',' -f2,3,4 \
    > ${1}.csv

# sed:
# replace spaces by ','
# stirp leading spaces
# remove first line of original file (nnode etc.)
# remove last line (comment from triangle)

# cut:
# get fields 2 and 3 (x and y)
