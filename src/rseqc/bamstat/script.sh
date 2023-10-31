#!/bin/bash

set -eo pipefail 

filename="$(basename -- $par_input/*.bam)"
mkdir -p $par_output

bam_stat.py --input $par_input/$filename > $par_output/${filename%%.*}.txt