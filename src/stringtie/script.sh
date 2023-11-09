#!/bin/bash

set -eo pipefail

if [ $par_strandedness=='forward' ]; then
    strand='--fr'
elif [ $par_strandedness=='reverse' ]; then
    strand='--rf'
fi

stringtie \
    $par_bam \
    $strand \
    ${par_annotation_gtf:+-G $par_annotation_gtf} \
    -o $par_transcript_gtf \
    -A "abundance.txt" \
    ${par_annotation_gtf:+-C "coverage.gtf"} \
    ${par_annotation_gtf:+-b "ballgown"} \
    -p $meta_cpus \
    $par_extra_stringtie_args \
    ${par_stringtie_ignore_gtf:+-e}

mv coverage.gtf $par_coverage_gtf
mv ballgown $par_ballgown
mv abundance.txt $par_abundance
