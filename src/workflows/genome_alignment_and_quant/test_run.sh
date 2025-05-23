#!/bin/bash

viash ns build --setup cb --parallel

# Split error message from standard output
# viash ns list > /dev/null 

echo "> Preparing reference data files"
gunzip --keep testData/minimal_test/reference/genes.gtf.gz
mkdir -p testData/minimal_test/reference/rsem_index
tar -C testData/minimal_test/reference/rsem_index --strip-components 1 -xavf testData/minimal_test/reference/rsem.tar.gz --no-same-owner

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
  --rsem_index testData/minimal_test/reference/rsem_index \
  --aligner star_rsem \
  -profile docker \
  -resume

echo "Removing reference data files"
rm testData/minimal_test/reference/genes.gtf
