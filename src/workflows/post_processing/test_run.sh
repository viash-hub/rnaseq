#!/bin/bash

viash ns build --setup cb --parallel

nextflow run target/nextflow/workflows/post_processing/main.nf \
  --publish_dir "testData/paired_end_test" \
  --id SRR6357070 \
  --paired true \
  --strandedness unstranded \
  --fasta "testData/test_output/reference_genome.fasta" \
  --fai "testData/test_output/reference_genome.fasta.fai" \
  --gtf "testData/test_output/gene_annotation.gtf" \
  --genome_bam "testData/paired_end_test/SRR6357070.genome.bam" \
  --chrom_sizes "testData/test_output/reference_genome.fasta.sizes" \
  --star_multiqc "testData/paired_end_test/SRR6357070.star_align.log" \
  --extra_picard_args "--ASSUME_SORTED true --REMOVE_DUPLICATES false --VALIDATION_STRINGENCY LENIENT --TMP_DIR tmp" \
  --gencode false \
  --biotype gene_biotype \
  -profile docker \
