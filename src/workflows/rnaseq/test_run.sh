#!/bin/bash

# viash ns build --setup cb --parallel

# Test single-end data
# cat > testData/minimal_test/sample_sheet.csv << HERE
# id,fastq_1
# SRR6357070_1,SRR6357070_1.fastq.gz
# HERE

# Test paired-end data
# cat > testData/minimal_test/sample_sheet.csv << HERE
# id,fastq_1,fastq_2,strandedness
# SRR6357070,SRR6357070_1.fastq.gz,SRR6357070_2.fastq.gz,reverse
# SRR6357071,SRR6357071_1.fastq.gz,SRR6357071_2.fastq.gz,reverse
# HERE

# nextflow run target/nextflow/workflows/rnaseq/main.nf \
#   --param_list testData/minimal_test/sample_sheet.csv \
#   --publish_dir "testData/full_pipeline_test" \
#   --umitools_bc_pattern "NNNN" \
#   --umitools_bc_pattern2 "NNNN" \
#   --fasta testData/reference/genome.fasta \
#   --gtf testData/reference/genes.gtf.gz \
#   --additional_fasta testData/reference/gfp.fa.gz \
#   --transcript_fasta testData/reference/transcriptome.fasta \
#   --bbsplit_fasta_list testData/reference/bbsplit_fasta_list.txt \
#   -profile docker \
#   -resume

  # nextflow run target/nextflow/workflows/rnaseq/main.nf \
  # --id SRR6357070 \
  # --fastq_1 testData/minimal_test/SRR6357070_1.fastq.gz \
  # --fastq_2 testData/minimal_test/SRR6357070_2.fastq.gz \
  # --publish_dir "testData/full_pipeline_test" \
  # --fasta testData/reference/genome.fasta \
  # --gtf testData/reference/genes.gtf.gz \
  # -profile docker 


cat > testData/minimal_test/sample_sheet.csv << HERE
id,fastq_1,fastq_2,strandedness
WT_REP1,SRR6357070_1.fastq.gz;SRR6357071_1.fastq.gz,SRR6357070_2.fastq.gz;SRR6357071_2.fastq.gz,reverse
WT_REP2,SRR6357072_1.fastq.gz,SRR6357072_2.fastq.gz,reverse
RAP1_IAA_30M_REP1,SRR6357076_1.fastq.gz,SRR6357076_2.fastq.gz,reverse
RAP1_UNINDUCED_REP1,SRR6357073_1.fastq.gz,,reverse
RAP1_UNINDUCED_REP2,SRR6357074_1.fastq.gz;SRR6357075_1.fastq.gz,,reverse
HERE

nextflow run target/nextflow/workflows/rnaseq/main.nf \
  --param_list testData/minimal_test/sample_sheet.csv \
  --publish_dir "testData/full_pipeline_test" \
  --fasta testData/reference/genome.fasta \
  --gtf testData/reference/genes.gtf.gz \
  --additional_fasta testData/reference/gfp.fa.gz \
  --transcript_fasta testData/reference/transcriptome.fasta \
  --bbsplit_fasta_list testData/reference/bbsplit_fasta_list.txt \
  --salmon_index testData/reference/salmon.tar.gz \
  -profile docker