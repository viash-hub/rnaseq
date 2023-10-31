#!/bin/bash

set -eo pipefail 

filename="$(basename -- $par_bam_input)"
ref_file="$(basename -- $par_bed_input)"
mkdir -p $par_output

infer_experiment.py -i $par_bam_input -r $par_bed_input > $par_output/${filename%%.*}.${ref_file%%.*}.infer_experiment.txt