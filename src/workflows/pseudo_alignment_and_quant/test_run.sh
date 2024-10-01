#!/bin/bash

# viash ns build --setup cb -q pseudo_alignment_and_quant

# Split error message from standard output
# viash ns list > /dev/null 

echo "> Preparing reference data files"
gunzip --keep testData/minimal_test/reference/genes.gtf.gz
mkdir -p testData/minimal_test/reference/salmon_index
tar -C testData/minimal_test/reference/salmon_index --strip-components 1 -xavf testData/minimal_test/reference/salmon.tar.gz 

cat > testData/minimal_test/input_fastq/sample_sheet.csv << HERE
id,fastq_1,fastq_2,strandedness
WT_REP1,SRR6357070_1.fastq.gz,SRR6357070_2.fastq.gz,reverse
RAP1_UNINDUCED_REP1,SRR6357073_1.fastq.gz,,reverse
HERE

echo "> Test 1: Salmon qunatification"
nextflow run target/nextflow/workflows/pseudo_alignment_and_quant/main.nf \
  --param_list testData/minimal_test/input_fastq/sample_sheet.csv \
  --publish_dir "test_results/pseudo_alignment_test1" \
  --fasta testData/minimal_test/reference/genome.fasta \
  --gtf testData/minimal_test/reference/genes.gtf.gz \
  --transcript_fasta testData/minimal_test/reference/transcriptome.fasta \
  --salmon_index testData/minimal_test/reference/salmon_index \
  --pseudo_aligner salmon \
  -profile docker \
  -resume

# echo "> Test 2: Kallisto qunatification"
# nextflow run target/nextflow/workflows/pseudo_alignment_and_quant/main.nf \
#   --param_list testData/minimal_test/input_fastq/sample_sheet.csv \
#   --publish_dir "test_results/pseudo_alignment_test2" \
#   --fasta testData/minimal_test/reference/genome.fasta \
#   --gtf testData/minimal_test/reference/genes.gtf.gz \
#   --transcript_fasta testData/minimal_test/reference/transcriptome.fasta \
#   --kallisto_index test_results/prepare_genome_test3/Kallisto_index \
#   --pseudo_aligner kallisto \
#   -profile docker \
#   -resume

echo "Removing reference data files"
rm testData/minimal_test/reference/genes.gtf
rm -r testData/minimal_test/reference/salmon_index
