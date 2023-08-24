#!/bin/bash

# bin/viash ns build --setup cb

nextflow run target/nextflow/workflows/pre_processing/main.nf \
  --id SRR6357070 \
  --input 'testData/test/SRR6357070_1.fastq.gz' \
  --publish_dir "testData/test_output" \
  -profile docker \
  -resume


