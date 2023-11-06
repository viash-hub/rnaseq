#!/bin/bash

set -eo pipefail 

mkdir -p $par_output

bam_stat.py \
    --input $par_input \
    --mapq $par_map_qual \
> $par_output