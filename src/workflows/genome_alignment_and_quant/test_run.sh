#!/bin/bash

viash ns build --setup cb --parallel

# Split error message from standard output
# viash ns list > /dev/null 

CURR=`pwd` 

cat > testData/minimal_test/input_fastq/sample_sheet.csv << HERE
id,fastq_1,fastq_2,strandedness
WT_REP1,SRR6357070_1.fastq.gz,SRR6357070_2.fastq.gz,reverse
WT_REP2,SRR6357072_1.fastq.gz,SRR6357072_2.fastq.gz,reverse
HERE

nextflow run target/nextflow/workflows/genome_alignment_and_quant/main.nf \
  --param_list testData/minimal_test/input_fastq/sample_sheet.csv \
  --publish_dir "test_results/genome_alignment_test" \
  --fasta testData/minimal_test/reference/genome.fasta \
  --gtf testData/minimal_test/reference/genes.gtf.gz \
  --transcript_fasta testData/minimal_test/reference/transcriptome.fasta \
  --star_index testData/test_output/star_index \
  --extra_star_align_args "--readFilesCommand gunzip -c --quantMode TranscriptomeSAM --twopassMode Basic --outSAMtype BAM Unsorted --runRNGseed 0 --outFilterMultimapNmax 20 --alignSJDBoverhangMin 1 --outSAMattributes NH HI AS NM MD --quantTranscriptomeBan Singleend --outSAMstrandField intronMotif" \
  -profile docker \
  # -resume