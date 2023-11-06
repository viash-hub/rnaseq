#!/bin/bash

set -exo pipefail 

filename="$(basename -- $par_input)"

if $par_paired; then
    inner_distance.py -i $par_input -r $par_refgene -o $filename > stdout.txt
    head -n 2 stdout.txt > $filename.inner_distance_stats.txt
fi
