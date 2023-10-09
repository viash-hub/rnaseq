#!/bin/bash

set -eo pipefail

## VIASH START
meta_cpus=2
## VIASH END

samtools sort -@ $meta_cpus -o $par_output $par_input/*d.out.bam