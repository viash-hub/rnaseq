#!/bin/bash

viash ns build --setup cb --parallel 

# echo "Test 1: Annotation file format - GTF"
# nextflow run target/nextflow/workflows/prepare_genome/main.nf \
#     --id test1 \
#     --publish_dir "test_results/prepare_genome_test1" \
#     --fasta testData/minimal_test/reference/genome.fasta \
#     --gtf testData/minimal_test/reference/genes.gtf.gz \
#     --additional_fasta testData/minimal_test/reference/gfp.fa.gz \
#     --transcript_fasta testData/minimal_test/reference/transcriptome.fasta \
#     --genotype false \
#     --biotype gene_biotype \
#     --bbsplit_fasta_list "testData/minimal_test/reference/bbsplit_fasta/sarscov2.fa;testData/minimal_test/reference/bbsplit_fasta/human.fa" \
#     --salmon_index testData/minimal_test/reference/salmon.tar.gz \
#     --rsem_index testData/minimal_test/reference/rsem.tar.gz \
#     -profile docker \
#     -resume

# echo "Test 2: Annotation file format - GFF"
# nextflow run target/nextflow/workflows/prepare_genome/main.nf \
#     --id test2 \
#     --publish_dir "test_results/prepare_genome_test2" \
#     --fasta testData/minimal_test/reference/genome.fasta \
#     --gff testData/minimal_test/reference/genes.gff.gz \
#     --additional_fasta testData/minimal_test/reference/gfp.fa.gz \
#     --transcript_fasta testData/minimal_test/reference/transcriptome.fasta \
#     --genotype false \
#     --biotype gene_biotype \
#     --bbsplit_fasta_list "testData/minimal_test/reference/bbsplit_fasta/sarscov2.fa;testData/minimal_test/reference/bbsplit_fasta/human.fa" \
#     --salmon_index testData/minimal_test/reference/salmon.tar.gz \
#     --rsem_index testData/minimal_test/reference/rsem.tar.gz \
#     -profile docker \
#     -resume

echo "Test 3: Annotation file format - GTF; Generate indices; Generate transcripts fasta"
nextflow run target/nextflow/workflows/prepare_genome/main.nf \
    --id test3 \
    --publish_dir "test_results/prepare_genome_test3" \
    --fasta testData/minimal_test/reference/genome.fasta \
    --gtf testData/minimal_test/reference/genes.gtf.gz \
    --additional_fasta testData/minimal_test/reference/gfp.fa.gz \
    --genotype false \
    --biotype gene_biotype \
    --bbsplit_fasta_list "testData/minimal_test/reference/bbsplit_fasta/sarscov2.fa;testData/minimal_test/reference/bbsplit_fasta/human.fa" \
    --pseudo_aligner kallisto \
    --aligner star_rsem \
    -profile docker \
    -resume
