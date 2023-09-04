#!/bin/bash

bin/viash ns build --setup cb

# nextflow run target/nextflow/workflows/pre_processing/main.nf \
#   --id SRR6357070_1 \
#   --input 'testData/test/SRR6357070_1.fastq.gz' \
#   --publish_dir "testData/test_output" \
#   --umitools_bc_pattern "NNNN" \
#   # -profile docker \
#   # -resume

# Test single-end data
# cat > testData/test/sample_sheet.csv << HERE
# id,fastq_1
# SRR6357070_1,SRR6357070_1.fastq.gz
# HERE

# nextflow run target/nextflow/workflows/pre_processing/main.nf \
#   --param_list testData/test/sample_sheet.csv \
#   --publish_dir "testData/single_end_test" \
#   --umitools_bc_pattern "NNNN"

# Test paired-end data
cat > testData/test/sample_sheet.csv << HERE
id,fastq_1,fastq_2
SRR6357070,SRR6357070_1.fastq.gz,SRR6357070_2.fastq.gz
HERE

nextflow run target/nextflow/workflows/pre_processing/main.nf \
  --param_list testData/test/sample_sheet.csv \
  --publish_dir "testData/paired_end_test" \
  --umitools_bc_pattern "NNNN" \
  --umitools_bc_pattern2 "NNNN" \
#   -profile docker \
#   -resume