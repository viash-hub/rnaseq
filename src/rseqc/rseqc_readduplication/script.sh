#!/bin/bash

set -eo pipefail 

mkdir -p $par_output

read_duplication.py -i $par_bam_input -o $par_output 
