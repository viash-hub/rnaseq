#!/bin/bash

# v;iash ns build --setup cb --parallel

# Split error message from standard output
# viash ns list > /dev/null 

echo "> Preparing reference data files"
gunzip --keep testData/minimal_test/reference/genes.gtf.gz

cat > testData/minimal_test/input_fastq/sample_sheet.csv << HERE
id,fastq_1,fastq_2,strandedness
WT_REP1,SRR6357070_1.fastq.gz,SRR6357070_2.fastq.gz,reverse
RAP1_UNINDUCED_REP1,SRR6357073_1.fastq.gz,,reverse
HERE

# echo "> Test 1: STAR Salmon"
# nextflow run target/nextflow/workflows/genome_alignment_and_quant/main.nf \
#   --param_list testData/minimal_test/input_fastq/sample_sheet.csv \
#   --publish_dir test_results/genome_alignment_test1 \
#   --fasta testData/minimal_test/reference/genome.fasta \
#   --gtf testData/minimal_test/reference/genes.gtf \
#   --transcript_fasta testData/minimal_test/reference/transcriptome.fasta \
#   --star_index test_results/output_test1/STAR_index \
#   --aligner star_salmon \
#   -profile docker \
#   -resume

echo "> Test 2: STAR RSEM"
nextflow run target/nextflow/workflows/genome_alignment_and_quant/main.nf \
  --param_list testData/minimal_test/input_fastq/sample_sheet.csv \
  --publish_dir test_results/genome_alignment_test2 \
  --fasta testData/minimal_test/reference/genome.fasta \
  --gtf testData/minimal_test/reference/genes.gtf \
  --transcript_fasta testData/minimal_test/reference/transcriptome.fasta \
  --rsem_index test_results/output_test1/RSEM_index \
  --aligner star_rsem \
  --extra_rsem_calculate_expression_args "--star --star-output-genome-bam --star-gzipped-read-file --estimate-rspd --seed 1" \
  -profile docker \
  -resume

echo "Removing reference data files"
rm testData/minimal_test/reference/genes.gtf
