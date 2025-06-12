nextflow run target/nextflow/prepare_reference/main.nf \
    --input_genome_fasta resources_test/minimal_test/reference/genome.fasta \
    --input_transcriptome_gtf resources_test/minimal_test/reference/genes.gtf.gz \
    --publish_dir test_results/prepare_genome \
    -profile docker