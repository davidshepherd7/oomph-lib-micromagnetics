#!/bin/bash

set -o errexit
set -o nounset

if [[ ! "$#" -eq 2 ]]; then
    echo "Usage: trace_headers.sh [tracefile] [header]"
    exit 0
fi

col_nums=$(cat "$1" | head -n1 | tr -s '; ' '\n' | grep -n "$2" --color=never \
		  | cut -d : -f 1 | head -c-1 | tr '\n' ',')

cut "$1" -d \; -f"$col_nums"
