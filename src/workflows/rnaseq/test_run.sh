#!/bin/bash

viash ns build --setup cb --parallel

# cat > testData/minimal_test/input_fastq/sample_sheet.csv << HERE
# id,fastq_1,fastq_2,strandedness
# WT_REP1,SRR6357070_1.fastq.gz;SRR6357071_1.fastq.gz,SRR6357070_2.fastq.gz;SRR6357071_2.fastq.gz,reverse
# WT_REP2,SRR6357072_1.fastq.gz,SRR6357072_2.fastq.gz,reverse
# RAP1_IAA_30M_REP1,SRR6357076_1.fastq.gz,SRR6357076_2.fastq.gz,reverse
# RAP1_UNINDUCED_REP1,SRR6357073_1.fastq.gz,,reverse
# RAP1_UNINDUCED_REP2,SRR6357074_1.fastq.gz;SRR6357075_1.fastq.gz,,reverse
# HERE

nextflow run target/nextflow/workflows/rnaseq/main.nf \
  --param_list testData/minimal_test/input_fastq/sample_sheet.csv \
  --publish_dir "test_results/full_pipeline_test" \
  --fasta testData/minimal_test/reference/genome.fasta \
  --gtf testData/minimal_test/reference/genes.gtf.gz \
  --additional_fasta testData/minimal_test/reference/gfp.fa.gz \
  --transcript_fasta testData/minimal_test/reference/transcriptome.fasta \
  --bbsplit_fasta_list testData/minimal_test/reference/bbsplit_fasta_list.txt \
  --salmon_index testData/minimal_test/reference/salmon.tar.gz \
  --skip_pseudo_alignment \
  -profile docker --resume