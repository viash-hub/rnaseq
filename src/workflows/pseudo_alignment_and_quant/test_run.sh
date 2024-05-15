#!/bin/bash

viash ns build --setup cb -q pseudo_alignment_and_quant

# Split error message from standard output
# viash ns list > /dev/null 

CURR=`pwd` 

# Test paired-end data
cat > testData/minimal_test/input_fastq/sample_sheet.csv << HERE
id,fastq_1,fastq_2,strandedness
WT_REP1,SRR6357070_1.fastq.gz,SRR6357070_2.fastq.gz,reverse
WT_REP2,SRR6357072_1.fastq.gz,SRR6357072_2.fastq.gz,reverse
HERE

nextflow run target/nextflow/workflows/pseudo_alignment_and_quant/main.nf \
  --param_list testData/minimal_test/input_fastq/sample_sheet.csv \
  --publish_dir "test_results/psudo_alignment_test" \
  --fasta testData/minimal_test/reference/genome.fasta \
  --gtf testData/minimal_test/reference/genes.gtf.gz \
  --transcript_fasta testData/minimal_test/reference/transcriptome.fasta \
  --salmon_index testData/minimal_test/reference/salmon.tar.gz \
  --pseudo_aligner salmon \
  -profile docker \
  # -resume