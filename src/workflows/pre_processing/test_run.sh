#!/bin/bash

viash ns build --parallel --setup cb

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
# SRR6357071_1,SRR6357071_1.fastq.gz
# HERE

# nextflow run target/nextflow/workflows/pre_processing/main.nf \
#   --param_list testData/test/sample_sheet.csv \
#   --publish_dir "testData/single_end_test" \
#   --umitools_bc_pattern "NNNN" \
#   --bbsplit_fasta_list testData/reference/bbsplit_fasta_list.txt \
#   --fasta testData/test_output/ref.prepare_genome.fasta_uncompressed \
#   --bbsplit_index testData/test_output/ref.prepare_genome.bbsplit_index_uncompressed \
#   --ribo_database_manifest testData/reference/rrna-db-defaults.txt 

# Test paired-end data
cat > testData/test/sample_sheet.csv << HERE
id,fastq_1,fastq_2
SRR6357070,SRR6357070_1.fastq.gz,SRR6357070_2.fastq.gz
SRR6357071,SRR6357071_1.fastq.gz,SRR6357071_2.fastq.gz
HERE

nextflow run target/nextflow/workflows/pre_processing/main.nf \
  --param_list testData/test/sample_sheet.csv \
  --publish_dir "testData/paired_end_test" \
  --umitools_bc_pattern "NNNN" \
  --umitools_bc_pattern2 "NNNN" \
  --bbsplit_fasta_list testData/reference/bbsplit_fasta_list.txt \
  --transcript_fasta testData/test_output/transcriptome.fasta \
  --gtf testData/test_output/gene_annotation.gtf \
  --bbsplit_index testData/test_output/bbsplit_index \
  --ribo_database_manifest testData/reference/rrna-db-defaults.txt \
  --salmon_index testData/test_output/salmon_index \
  --skip_trimming false \
  --remove_ribo_rna true
  # -profile docker \
#   -resume