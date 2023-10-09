#!/bin/bash

set -eo pipefail

## VIASH START
meta_cpus=2
## VIASH END

if [ -f "$par_fasta" ]; then
    reference="--reference $par_fasta"
else
    reference=""
fi

samtools stats --threads $meta_cpus $reference $par_input/*.bam > $par_output