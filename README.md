# RNAseq.vsh

<!-- README.md is generated by running 'quarto render README.qmd' -->

A version of the [nf-core/rnaseq](https://github.com/nf-core/rnaseq)
pipeline (version 3.12.0) in the [Viash framework](http://www.viash.io).

## Rationale

We stick to the original nf-core pipeline as much as possible. This also
means that we create a subworkflow for the 5 main stages of the pipeline
as depicted in the [README](https://github.com/nf-core/rnaseq).

The original version of the nf-core/rnaseq pipeline allowed one to point
to input files directly, but more recently [a sample sheet file is
required](https://github.com/nf-core/rnaseq#usage). We don’t have to
implement the sample sheet approach just yet and just take fastq files
as input for the moment.

## Getting started

As test data, we can use the small dataset nf-core provided with [their
`test`
profile](https://github.com/nf-core/test-datasets/blob/rnaseq3/samplesheet/v3.10/samplesheet_test.csv):
<https://github.com/nf-core/test-datasets/tree/rnaseq3/testdata/GSE110004>.

A simple script has been provided to fetch those files from the github
repository and store them under `testData/test` (the subdirectory is
created to support `full_test` later as well): `bin/get_testData.sh`.
Additional test data are fetched from (provided by rseqc package) to
evaluate paired-end experiments.

Additionally, a script has been provided to fetch all the necessary
reference data files from the github repository
(https://github.com/nf-core/test-datasets/tree/rnaseq3/reference) and
store them under `testData/reference`: `bin/get_reference.sh`

To get started, we need to:

1.  [Install
    `nextflow`](https://www.nextflow.io/docs/latest/getstarted.html)
    system-wide

2.  Fetch the test data:

``` bash
bin/get_testData.sh
bin/get_reference.sh
```

## Running the pipeline

To actually run the pipeline, we first need to build the components and
pipeline:

``` bash
viash ns build --setup cb
```

Now we can run the pipeline using the command:

``` bash
nextflow run target/nextflow/workflows/pre_processing/main.nf \
  -profile docker \
  --id test \
  --input testData/test/SRR6357070_1.fastq.gz \
  --publish_dir testData/test_output/
```

Alternatively, we can run the pipeline with a sample sheet using the
built-in `--param_list` functionality: (Read file paths must be
specified relative to the sample sheet’s path)

``` bash
cat > testData/test/sample_sheet.csv << HERE
id,fastq_1
SRR6357070_1,SRR6357070_1.fastq.gz
SRR6357071_1,SRR6357071_1.fastq.gz
HERE

nextflow run target/nextflow/workflows/pre_processing/main.nf \
  --param_list testData/test/sample_sheet.csv \
  --publish_dir "testData/test_output" \
  -profile docker \
  -resume
```

## Pipeline sub-workflows and components

The pipeline has 5 sub-workflows that can be run separately.

1.  Prepare genome: This sub-workflow can be used to prepare all the
    necessary reference data, i.e., uncompress provided reference data
    or generate required index files.

The following command can be used to run the workflow:

``` bash
nextflow run target/nextflow/workflows/prepare_genome/main.nf \
    --id id \
    --publish_dir "testData/test_output" \
    --fasta testData/reference/genome.fasta \
    --gtf testData/reference/genes.gtf.gz \
    --additional_fasta testData/reference/gfp.fa.gz \
    --transcript_fasta testData/reference/transcriptome.fasta \
    --gencode true \
    --biotype gene_type \
    --bbsplit_fasta_list testData/reference/bbsplit_fasta_list.txt
```

2.  Pre-processing: This sub-workflow can be used perform QC on the
    RNA-seq data. At the moment, it performs FastQC, extracts UMIs,
    trims adapters, and removes ribosomal RNA reads.

The following command can be used to run the workflow:

``` bash
cat > testData/test/sample_sheet.csv << HERE
id,fastq_1,fastq_2
SRR6357070,SRR6357070_1.fastq.gz,SRR6357070_2.fastq.gz
HERE

nextflow run target/nextflow/workflows/pre_processing/main.nf \
  --param_list testData/test/sample_sheet.csv \
  --publish_dir "testData/paired_end_test" \
  --umitools_bc_pattern "NNNN" \
  --umitools_bc_pattern2 "NNNN" \
  --bbsplit_fasta_list testData/reference/bbsplit_fasta_list.txt \
  --fasta testData/test_output/ref.prepare_genome.fasta_uncompressed \
  --bbsplit_index testData/test_output/ref.prepare_genome.bbsplit_index_uncompressed \
  --ribo_database_manifest testData/reference/rrna-db-defaults.txt
```

3.  Genome alignment and quantification: This sub-workflow can be used
    to perform genome alignment using STAR and alignment quantification
    using Salmon. Alignment sorting and indexing, as well as as well as
    computation of statistics from the BAM files is performed using
    Samtools. UMI-based deduplication can also be performed.

The following command can be used to run the workflow:

``` bash
cat > testData/test/sample_sheet.csv << HERE
id,fastq_1,fastq_2
SRR6357070,SRR6357070_1.fastq.gz,SRR6357070_2.fastq.gz
HERE

nextflow run target/nextflow/workflows/genome_alignment_and_quant/main.nf \
  --param_list testData/test/sample_sheet.csv \
  --publish_dir "testData/paired_end_test" \
  --fasta testData/test_output/ref.prepare_genome.fasta_uncompressed \
  --gtf testData/test_output/ref.prepare_genome.gtf_uncompressed.gtf \
  --star_index testData/test_output/ref.prepare_genome.star_index_uncompressed \
  --transcript_fasta testData/test_output/ref.prepare_genome.transcript_fasta_uncompressed.fasta \
  --extra_star_align_args "--readFilesCommand gunzip -c --quantMode TranscriptomeSAM --twopassMode Basic --outSAMtype BAM Unsorted --runRNGseed 0 --outFilterMultimapNmax 20 --alignSJDBoverhangMin 1 --outSAMattributes NH HI AS NM MD --quantTranscriptomeBan Singleend --outSAMstrandField intronMotif"
```

4.  Post-processing: This sub-workflow can be used for duplicate read
    marking (picard MarkDuplicates), transcript assembly and
    quantification (StringTie), feature biotype QC (featureCounts) and
    creation of bigWig coverage files.

5.  Final QC: This sub-workflow performs extensive quality control
    (RSeQC, dupRadar, Qualimap, Preseq, DESeq2) and presents QC for raw
    read, alignment, gene biotype, sample similarity, and strand
    specificity (MultiQC).
