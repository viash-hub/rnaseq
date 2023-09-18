#!/bin/bash

bin/viash ns build --setup cb

# Test single-end data
cat > testData/test/sample_sheet.csv << HERE
id,fastq_1
SRR6357070_1,SRR6357070_1.fastq.gz
HERE

nextflow run target/nextflow/workflows/rnaseq/main.nf \
  --param_list testData/test/sample_sheet.csv \
  --publish_dir "testData/single_end_test" \
  --umitools_bc_pattern "NNNN" \
  --fasta testData/reference/genome.fasta \
  --gtf testData/reference/genes.gtf.gz \
  --additional_fasta testData/reference/gfp.fa.gz \
  --transcript_fasta testData/reference/transcriptome.fasta \
  --gencode
  # -profile docker \
  # -resume

  # --gff testData/reference/genes.gff.gz \
  # --additional_fasta testData/reference/gfp.fa.gz \
  # --transcript_fasta testData/reference/transcriptome.fasta \
  # --bbsplit_fasta_list testData/reference/bbsplit_fasta_list.txt \
  # --rsem_index testData/reference/rsem.tar.gz \
  # --salmon_index testData/reference/salmon.tar.gz \
  # --hisat2_index testData/reference/hisat2.tar.gz \
  # --gencode true \
  # --biotype gene_type \

# Test paired-end data
# cat > testData/test/sample_sheet.csv << HERE
# id,fastq_1,fastq_2
# SRR6357070,SRR6357070_1.fastq.gz,SRR6357070_2.fastq.gz
# HERE

# nextflow run target/nextflow/workflows/pre_processing/main.nf \
#   --param_list testData/test/sample_sheet.csv \
#   --publish_dir "testData/paired_end_test" \
#   --umitools_bc_pattern "NNNN" \
#   --umitools_bc_pattern2 "NNNN" \
#   -profile docker \
#   -resume