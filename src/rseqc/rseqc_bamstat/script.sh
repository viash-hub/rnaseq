#!/bin/bash

set -eo pipefail 

filename="$(basename -- $par_bam_input)"
mkdir -p $par_output

bam_stat.py --input $par_bam_input > $par_output/${filename%%.*}.bamstat.txt