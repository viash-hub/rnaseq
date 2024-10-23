#!/bin/bash

viash ns build --parallel --setup cb #-q pre_processing

echo "> Preparing reference data files"
gunzip --keep testData/minimal_test/reference/genes.gtf.gz
mkdir -p testData/minimal_test/reference/salmon_index
tar -C testData/minimal_test/reference/salmon_index --strip-components 1 -xavf testData/minimal_test/reference/salmon.tar.gz --no-same-owner

# Test paired-end data
cat > testData/minimal_test/input_fastq/sample_sheet.csv << HERE
id,fastq_1,fastq_2,strandedness
WT_REP2,SRR6357072_1.fastq.gz,SRR6357072_2.fastq.gz,auto
RAP1_UNINDUCED_REP1,SRR6357073_1.fastq.gz,,reverse
HERE

echo "> Test 1: Running workflow with trimgalore"
nextflow run target/nextflow/workflows/pre_processing/main.nf \
  --param_list testData/minimal_test/input_fastq/sample_sheet.csv \
  --publish_dir "test_results/pre_processing_test1" \
  --bbsplit_fasta_list testData/minimal_test/reference/bbsplit_fasta_list.txt \
  --transcript_fasta testData/minimal_test/reference/transcriptome.fasta \
  --gtf testData/minimal_test/reference/genes.gtf \
  --salmon_index testData/minimal_test/reference/salmon_index \
  --skip_trimming false \
  --trimmer trimgalore \
  --remove_ribo_rna true \
  --ribo_database_manifest testData/minimal_test/reference/rrna-db-defaults.txt \
  --skip_bbsplit true \
  --bbsplit_index test_results/prepare_genome_test1/BBSplit_index \
  --with_umi false \
  -profile docker \
  -resume

# echo "> Test 2: Running workflow with fastp"
# nextflow run target/nextflow/workflows/pre_processing/main.nf \
#   --param_list testData/minimal_test/input_fastq/sample_sheet.csv \
#   --publish_dir "test_results/pre_processing_test2" \
#   --bbsplit_fasta_list testData/minimal_test/reference/bbsplit_fasta_list.txt \
#   --transcript_fasta testData/minimal_test/reference/transcriptome.fasta \
#   --gtf testData/minimal_test/reference/genes.gtf \
#   --salmon_index testData/minimal_test/reference/salmon_index \
#   --skip_trimming false \
#   --trimmer fastp \
#   --remove_ribo_rna false \
#   --ribo_database_manifest testData/minimal_test/reference/rrna-db-defaults.txt \
#   --skip_bbsplit false \
#   --bbsplit_index test_results/output_test1/BBSplit_index \
#   -profile docker \
#   -resume

echo "Removing reference data files"
rm testData/minimal_test/reference/genes.gtf
rm -r testData/minimal_test/reference/salmon_index

# TODO: Fix error while running sortmerna component
# docker: Error response from daemon: failed to create task for container: failed to create shim task: OCI runtime create failed: runc create failed: unable to start container process: exec: "/bin/bash": stat /bin/bash: no such file or directory: unknown.