#!/bin/bash

set -eo pipefail 
set -eo pipefail 

filename="$(basename -- $par_bam_input)"
ref_file="$(basename -- $par_bed_input)"
mkdir -p $par_output

junction_annotation.py -i $par_bam_input -r $par_bed_input -o $par_output 2> $par_output/${filename%%.*}.${ref_file%%.*}.junction_annotation.log
