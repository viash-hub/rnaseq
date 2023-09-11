#!/bin/bash

## VIASH START
meta_cpus=8
# meta_memory_b='123'
## VIASH END

set -eo pipefail

# args=if RSEM_PREPAREREFERENCE_GENOME then '--STAR' else if MAKE_TRANSCRIPTS_FASTA then ''
if [ $par_star ]; then 
    # memory = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    STAR --runMode genomeGenerate --genomeDir $par_rsem \
    --genomeFastaFiles $par_fasta --sjdbGTFfile $par_gtf \
    --runThreadN $meta_cpus # $meta_memory 
    
    rsem-prepare-reference --gtf $par_gtf --star \
    --num-threads $meta_cpus \
    $par_fasta "$par_rsem/genome"
    
    cp "$par_rsem/genome.transcripts.fa" $par_transcript_fasta

else 
    rsem-prepare-reference --gtf $par_gtf \
    --num-threads $meta_cpus \
    $par_fasta "$par_rsem/genome"
    
    cp "$par_rsem/genome.transcripts.fa" $par_transcript_fasta
fi

        