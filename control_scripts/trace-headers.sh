#!/bin/bash

set -o errexit
set -o nounset

if [[ ! "$#" -eq 2 ]]; then
    echo "Usage: trace_headers.sh [tracefile] [header]"
    exit 0
fi

cat "$1" | head -n1 | tr -s '; ' '\n' | grep -n "$2" --color=always
