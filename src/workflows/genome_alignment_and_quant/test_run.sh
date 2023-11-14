#!/bin/bash

viash ns build --setup cb --parallel

# Split error message from standard output
# viash ns list > /dev/null 

CURR=`pwd` 

# Test single-end data
# cat > testData/test/sample_sheet.csv << HERE
# id,fastq_1
# SRR6357070_1,SRR6357070_1.fastq.gz
# HERE

# nextflow run target/nextflow/workflows/genome_alignment_and_quant/main.nf \
#   --param_list testData/test/sample_sheet.csv \
#   --publish_dir "testData/single_end_test" \
#   --fasta testData/test_output/ref.prepare_genome.fasta_uncompressed \
#   --gtf testData/test_output/ref.prepare_genome.gtf_uncompressed.gtf \
#   --star_index testData/test_output/ref.prepare_genome.star_index_uncompressed \
#   --transcript_fasta testData/test_output/ref.prepare_genome.transcript_fasta_uncompressed.fasta \
#   --extra_star_align_args "--readFilesCommand gunzip -c --quantMode TranscriptomeSAM --twopassMode Basic --outSAMtype BAM Unsorted --runRNGseed 0 --outFilterMultimapNmax 20 --alignSJDBoverhangMin 1 --outSAMattributes NH HI AS NM MD --quantTranscriptomeBan Singleend --outSAMstrandField intronMotif" \
  # -profile docker \
  # -resume

# Test paired-end data
cat > testData/test/sample_sheet.csv << HERE
id,fastq_1,fastq_2
SRR6357070,$CURR/testData/paired_end_test/SRR6357070.pre_processing.read1.fq.gz,$CURR/testData/paired_end_test/SRR6357070.pre_processing.read2.fq.gz
SRR6357071,$CURR/testData/paired_end_test/SRR6357071.pre_processing.read1.fq.gz,$CURR/testData/paired_end_test/SRR6357071.pre_processing.read2.fq.gz
HERE

nextflow run target/nextflow/workflows/genome_alignment_and_quant/main.nf \
  --param_list testData/test/sample_sheet.csv \
  --publish_dir "testData/paired_end_test" \
  --fasta testData/test_output/reference_genome.fasta \
  --gtf testData/test_output/gene_annotation.gtf \
  --star_index testData/test_output/star_index \
  --transcript_fasta testData/test_output/transcriptome.fasta \
  --extra_star_align_args "--readFilesCommand gunzip -c --quantMode TranscriptomeSAM --twopassMode Basic --outSAMtype BAM Unsorted --runRNGseed 0 --outFilterMultimapNmax 20 --alignSJDBoverhangMin 1 --outSAMattributes NH HI AS NM MD --quantTranscriptomeBan Singleend --outSAMstrandField intronMotif" \
  # -profile docker \
  # -resume