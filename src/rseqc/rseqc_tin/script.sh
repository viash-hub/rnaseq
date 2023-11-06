#!/bin/bash

set -eo pipefail 

filename="$(basename -- $par_input)"

tin.py \
    -i $par_input \
    -r $par_refgene \
    -c $par_minimum_coverage \

mv ${filename%.*}.summary.txt $par_output_tin_summary
mv ${filename%.*}.tin.xls $par_output_tin
