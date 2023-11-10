#!/bin/bash

set -eo pipefail 

filename="$(basename -- $par_input)"

tin.py \
    -i $par_input \
    -r $par_refgene \
    -c $par_minimum_coverage \
    -n $par_sample_size \
    -s $par_subtract_background

mv ${filename%.*}.summary.txt $par_output_tin_summary
mv ${filename%.*}.tin.xls $par_output_tin
