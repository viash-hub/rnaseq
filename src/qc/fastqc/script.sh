#!/bin/bash

set -eo pipefail

mkdir -p $par_output 

if [ -f "$par_input_r2" ]; then
    # paired-end data
    eval fastqc $par_input $par_input_r2 -o $par_output
else
    # single-end data
    eval fastqc $par_input -o $par_output
fi


