#!/bin/bash

set -eo pipefail 

filename="$(basename -- $par_bam_input)"
ref_file="$(basename -- $par_refgene)"
mkdir -p $par_output

infer_experiment.py -i $par_bam_input -r > $par_output/${filename%%.*}.${ref_file%%.*}.infer_experiment.txt