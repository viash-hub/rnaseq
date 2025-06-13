#!/bin/bash

export NXF_VER=24.04.5

viash ns build

nextflow run . \
  -main-script target/nextflow/prepare_reads/main.nf \
  --input_r1 resources_test/minimal_test/input_fastq/SRR6357070_1.fastq.gz \
  --input_r2 resources_test/minimal_test/input_fastq/SRR6357070_2.fastq.gz \
  --publish_dir test_results/test_prepare_reads \
  -profile docker
