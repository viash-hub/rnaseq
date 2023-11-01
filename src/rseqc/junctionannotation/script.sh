#!/bin/bash

set -eo pipefail 
set -eo pipefail 

filename="$(basename -- $par_bam_input)"
ref_file="$(basename -- $par_refgene)"
mkdir -p $par_output

junction_annotation.py -i $par_bam_input -r $par_refgene -o $par_output 2> $par_output/${filename%%.*}.${ref_file%%.*}.junction_annotation.log
