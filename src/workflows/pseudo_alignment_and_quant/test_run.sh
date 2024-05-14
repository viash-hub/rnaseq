#!/bin/bash

viash ns build --setup cb --parallel

# Split error message from standard output
# viash ns list > /dev/null 

CURR=`pwd` 

# Test paired-end data
cat > testData/minimal_test/sample_sheet.csv << HERE
id,fastq_1,fastq_2
SRR6357070,$CURR/testData/paired_end_test/SRR6357070.pre_processing.read_1.fq.gz,$CURR/testData/paired_end_test/SRR6357070.pre_processing.read_2.fq.gz
SRR6357071,$CURR/testData/paired_end_test/SRR6357071.pre_processing.read_1.fq.gz,$CURR/testData/paired_end_test/SRR6357071.pre_processing.read_2.fq.gz
HERE

nextflow run target/nextflow/workflows/genome_alignment_and_quant/main.nf \
  --param_list testData/minimal_test/sample_sheet.csv \
  --publish_dir "testData/paired_end_test" \
  --fasta testData/test_output/reference_genome.fasta \
  --gtf testData/test_output/ \
  --transcript_fasta testData/test_output/transcriptome.fasta \
  --salmon_index testData/minimal_test/reference/salmon.tar.gz \
  --pseudo_aligner salmon
  # -profile docker \
  # -resume