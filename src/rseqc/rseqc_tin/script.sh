#!/bin/bash

set -eo pipefail 

filename="$(basename -- $par_bam_input)"
ref_file="$(basename -- $par_refgene)"
mkdir -p $par_output

read_distribution.py -i $par_bam_input -r $par_refgene -o $par_output