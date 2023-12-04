#!/bin/bash

bin/viash ns build --setup cb --parallel

# Test single-end data
# cat > testData/minimal_test/sample_sheet.csv << HERE
# id,fastq_1
# SRR6357070_1,SRR6357070_1.fastq.gz
# HERE

# Test paired-end data
cat > testData/minimal_test/sample_sheet.csv << HERE
id,fastq_1,fastq_2,strandedness
SRR6357070,SRR6357070_1.fastq.gz,SRR6357070_2.fastq.gz,reverse
SRR6357071,SRR6357071_1.fastq.gz,SRR6357071_2.fastq.gz,reverse
HERE

nextflow run target/nextflow/workflows/rnaseq/main.nf \
  --param_list testData/minimal_test/sample_sheet.csv \
  --publish_dir "testData/full_pipeline_test" \
  --umitools_bc_pattern "NNNN" \
  --umitools_bc_pattern2 "NNNN" \
  --fasta testData/reference/genome.fasta \
  --gtf testData/reference/genes.gtf.gz \
  --additional_fasta testData/reference/gfp.fa.gz \
  --transcript_fasta testData/reference/transcriptome.fasta \
  --additional_fasta testData/reference/gfp.fa.gz \
  --bbsplit_fasta_list testData/reference/bbsplit_fasta_list.txt \
  -profile docker \
  -resume