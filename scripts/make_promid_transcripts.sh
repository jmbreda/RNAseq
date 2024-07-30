#!/bin/bash

promoterome=$1
outfile=$2

cut -f9 "$promoterome" | cut -f1,3 -d';' | sed 's/ID=//g' | sed 's/;Members=/\t/g' | sed 's/,/\n\t/g' > "$outfile"