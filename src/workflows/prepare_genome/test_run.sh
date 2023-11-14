#!/bin/bash

viash ns build --setup cb --parallel

nextflow run target/nextflow/workflows/prepare_genome/main.nf \
    --id ref \
    --publish_dir "testData/test_output" \
    --fasta testData/reference/genome.fasta \
    --gtf testData/reference/genes.gtf.gz \
    --additional_fasta testData/reference/gfp.fa.gz \
    --transcript_fasta testData/reference/transcriptome.fasta \
    --genotype false \
    --biotype gene_biotype \
    --bbsplit_fasta_list testData/reference/bbsplit_fasta_list.txt \
    --salmon_index testData/reference/salmon.tar.gz \
    # -profile docker \ 
    # -resume
    # --gff testData/reference/genes.gff.gz \
    # --prepare_tools_indices a,b,c \
    # --gene_bed "" \
    # --splicesites "" \
    # --star_index "" \
    # --bbsplit_index "" \
    # --rsem_index testData/reference/rsem.tar.gz \
    # --salmon_index testData/reference/salmon.tar.gz \
    # --hisat2_index testData/reference/hisat2.tar.gz \
