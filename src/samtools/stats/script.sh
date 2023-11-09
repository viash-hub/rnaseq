#!/bin/bash

set -eo pipefail

if [ -f "$par_fasta" ]; then
    reference="--reference $par_fasta"
else
    reference=""
fi

samtools stats --threads $meta_cpus $reference $par_bam > $par_output