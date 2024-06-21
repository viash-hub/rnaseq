#!/bin/bash

set -eo pipefail

if $par_star; then 
    STAR --runMode genomeGenerate \
        --genomeDir $par_rsem \
        --genomeFastaFiles $par_fasta \
        --sjdbGTFfile $par_gtf \
        ${meta_cpus:+--runThreadN $meta_cpus} \
        ${meta_memory:+--limitGenomeGenerateRAM $meta_memory}
    rsem-prepare-reference \
        --gtf $par_gtf \
        --star \
        ${meta_cpus:+--num-threads $meta_cpus} \
        $par_fasta \
        $par_rsem/genome
    cp "$par_rsem/genome.transcripts.fa" $par_transcript_fasta
else 
    rsem-prepare-reference \
    --gtf $par_gtf \
    ${meta_cpus:+--num-threads $meta_cpus} \
    $par_fasta \
    "$par_rsem/genome"
    cp "$par_rsem/genome.transcripts.fa" $par_transcript_fasta
fi
    