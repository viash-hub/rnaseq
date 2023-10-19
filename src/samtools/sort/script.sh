#!/bin/bash

set -eo pipefail

## VIASH START
meta_cpus=2
## VIASH END

mkdir -p $par_output
filename="$(basename -- $par_input/*.bam)"
samtools sort -@ $meta_cpus -o $par_output/sorted_$filename $par_input/$filename