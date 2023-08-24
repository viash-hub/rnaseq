#!/bin/bash

# bin/viash ns build --setup cb

nextflow run target/nextflow/workflows/pre_processing/main.nf \
  --id SRR6357070_1 \
  --input 'testData/test/SRR6357070_1.fastq.gz' \
  --publish_dir "testData/test_output" \
  -profile docker \
  -resume


# cat > testData/test/sample_sheet.csv << HERE
# id,input
# SRR6357070_1,SRR6357070_1.fastq.gz
# SRR6357071_1,SRR6357071_1.fastq.gz
# HERE

# nextflow run target/nextflow/workflows/pre_processing/main.nf \
#   --param_list testData/test/sample_sheet.csv \
#   --publish_dir "testData/test_output" \
#   -profile docker \
#   -resume