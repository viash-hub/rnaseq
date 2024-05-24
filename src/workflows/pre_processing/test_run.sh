#!/bin/bash

viash ns build --parallel --setup cb 

# nextflow run target/nextflow/workflows/pre_processing/main.nf \
#   --id RAP1_UNINDUCED_REP1 \
#   --input 'testData/minimal_test/input_fastq/SRR6357073_1.fastq.gz' \
#   --publish_dir "test_results/preprocessing_no_samplesheet" \
#   --umitools_bc_pattern "NNNN" \
#   -profile docker \
#   -resume

# Test paired-end data
cat > testData/minimal_test/input_fastq/sample_sheet.csv << HERE
id,fastq_1,fastq_2,strandedness
WT_REP2,SRR6357072_1.fastq.gz,SRR6357072_2.fastq.gz,reverse
RAP1_IAA_30M_REP1,SRR6357076_1.fastq.gz,SRR6357076_2.fastq.gz,reverse
HERE

nextflow run target/nextflow/workflows/pre_processing/main.nf \
  --param_list testData/minimal_test/input_fastq/sample_sheet.csv \
  --publish_dir "testData/paired_end_test" \
  --umitools_bc_pattern "NNNN" \
  --umitools_bc_pattern2 "NNNN" \
  --bbsplit_fasta_list testData/minimal_test/reference/bbsplit_fasta_list.txt \
  --transcript_fasta testData/minimal_test/reference/transcriptome.fasta \
  --gtf testData/minimal_test/reference/gene_annotation.gtf \
  --salmon_index testData/minimal_test/reference/salmon.tar.gz \
  --skip_trimming false \
  --trimmer fastp \
  --remove_ribo_rna false \
  -profile docker \
  -resume