#!/bin/bash

set -eo pipefail 

filename="$(basename -- $par_bam_input)"
ref_file="$(basename -- $par_refgene)"
mkdir -p $par_output

if $par_paired; then
    inner_distance.py -i $par_bam_input -r $par_refgene -o $par_output > stdout.txt
    head -n 2 stdout.txt > $par_output/${filename%%.*}.${ref_file%%.*}.innerdistance.txt
fi
