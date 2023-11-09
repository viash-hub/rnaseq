#!/bin/bash

set -eo pipefail

avail_mem=3072
if [ ! $meta_memory_mb ]; then
    echo '[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
else
    avail_mem=$(( $meta_memory_mb*0.8 ))
fi

java -Xmx${avail_mem}M -jar $PICARD MarkDuplicates \
    $par_extra_picard_args \
    --INPUT $par_bam \
    --OUTPUT $par_output_bam \
    --REFERENCE_SEQUENCE $par_fasta \
    --METRICS_FILE $par_metrics