#!/bin/bash

set -eo pipefail

## VIASH START
par_input="test"
par_output="test_sorted"
meta_cpus=2
## VIASH END

samtools sort -@ $meta_cpus -o "$par_output.bam" -T $par_output $par_input/*d.out.bam